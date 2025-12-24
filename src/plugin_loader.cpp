/**
 * Plugin Loader Implementation
 *
 * Handles dynamic loading of VEP annotation plugins using dlopen/dlsym.
 */

#include "plugin.hpp"
#include "vep_annotator.hpp"
#include <dlfcn.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sstream>
#include <algorithm>

namespace vep {

PluginLoader::~PluginLoader() {
    // Unload all plugins in reverse order
    for (auto it = plugins_.rbegin(); it != plugins_.rend(); ++it) {
        unload_plugin(*it);
    }
    plugins_.clear();
}

bool PluginLoader::load_plugin(const std::string& path, const std::string& config) {
    // Clear any previous error
    dlerror();

    // Open the shared library
    void* handle = dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
    if (!handle) {
        last_error_ = "Failed to load plugin: " + std::string(dlerror());
        log(LogLevel::ERROR, last_error_);
        return false;
    }

    // Find the create_plugin function
    CreatePluginFunc create_func = reinterpret_cast<CreatePluginFunc>(
        dlsym(handle, "create_plugin")
    );
    if (!create_func) {
        last_error_ = "Plugin missing create_plugin function: " + path;
        log(LogLevel::ERROR, last_error_);
        dlclose(handle);
        return false;
    }

    // Find the destroy_plugin function
    DestroyPluginFunc destroy_func = reinterpret_cast<DestroyPluginFunc>(
        dlsym(handle, "destroy_plugin")
    );
    if (!destroy_func) {
        last_error_ = "Plugin missing destroy_plugin function: " + path;
        log(LogLevel::ERROR, last_error_);
        dlclose(handle);
        return false;
    }

    // Create the plugin instance
    VEPPlugin* plugin = nullptr;
    try {
        plugin = create_func();
    } catch (const std::exception& e) {
        last_error_ = "Plugin creation failed: " + std::string(e.what());
        log(LogLevel::ERROR, last_error_);
        dlclose(handle);
        return false;
    }

    if (!plugin) {
        last_error_ = "Plugin creation returned null: " + path;
        log(LogLevel::ERROR, last_error_);
        dlclose(handle);
        return false;
    }

    // Initialize the plugin
    PluginConfig plugin_config;
    plugin_config.plugin_path = path;
    plugin_config.config_string = config;

    // Parse config string (format: key1=value1;key2=value2)
    std::istringstream iss(config);
    std::string pair;
    while (std::getline(iss, pair, ';')) {
        size_t eq = pair.find('=');
        if (eq != std::string::npos) {
            std::string key = pair.substr(0, eq);
            std::string value = pair.substr(eq + 1);
            plugin_config.options[key] = value;
        }
    }

    // Validate configuration
    std::string validation_error = plugin->validate_config(plugin_config);
    if (!validation_error.empty()) {
        last_error_ = "Plugin configuration invalid: " + validation_error;
        log(LogLevel::ERROR, last_error_);
        destroy_func(plugin);
        dlclose(handle);
        return false;
    }

    // Initialize
    if (!plugin->initialize(plugin_config)) {
        last_error_ = "Plugin initialization failed: " + path;
        log(LogLevel::ERROR, last_error_);
        destroy_func(plugin);
        dlclose(handle);
        return false;
    }

    // Store the loaded plugin
    LoadedPlugin lp;
    lp.handle = handle;
    lp.plugin = plugin;
    lp.destroy_func = destroy_func;
    lp.path = path;
    plugins_.push_back(lp);

    auto info = plugin->get_info();
    log(LogLevel::INFO, "Loaded plugin: " + info.name + " v" + info.version);

    return true;
}

void PluginLoader::add_plugin_dir(const std::string& dir) {
    plugin_dirs_.push_back(dir);
}

void PluginLoader::load_plugins_from_dirs() {
    for (const auto& dir : plugin_dirs_) {
        DIR* d = opendir(dir.c_str());
        if (!d) continue;

        struct dirent* entry;
        while ((entry = readdir(d)) != nullptr) {
            std::string name = entry->d_name;

            // Check for shared library extension
#ifdef __APPLE__
            if (name.length() > 6 && name.substr(name.length() - 6) == ".dylib") {
#else
            if (name.length() > 3 && name.substr(name.length() - 3) == ".so") {
#endif
                std::string path = dir + "/" + name;

                // Check if it's a regular file
                struct stat st;
                if (stat(path.c_str(), &st) == 0 && S_ISREG(st.st_mode)) {
                    load_plugin(path);
                }
            }
        }

        closedir(d);
    }
}

std::vector<VEPPlugin*> PluginLoader::get_plugins() const {
    std::vector<VEPPlugin*> result;
    for (const auto& lp : plugins_) {
        result.push_back(lp.plugin);
    }
    return result;
}

std::vector<std::shared_ptr<AnnotationSource>> PluginLoader::get_all_sources() {
    std::vector<std::shared_ptr<AnnotationSource>> result;

    for (auto& lp : plugins_) {
        try {
            auto sources = lp.plugin->create_sources();
            result.insert(result.end(), sources.begin(), sources.end());
        } catch (const std::exception& e) {
            log(LogLevel::ERROR, "Plugin source creation failed: " + std::string(e.what()));
        }
    }

    return result;
}

void PluginLoader::unload_plugin(LoadedPlugin& p) {
    if (p.plugin && p.destroy_func) {
        try {
            p.destroy_func(p.plugin);
        } catch (...) {
            log(LogLevel::WARNING, "Plugin destruction failed: " + p.path);
        }
    }
    if (p.handle) {
        dlclose(p.handle);
    }
    p.plugin = nullptr;
    p.handle = nullptr;
}

} // namespace vep
