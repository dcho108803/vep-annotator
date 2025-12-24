/**
 * VEP Plugin Interface
 *
 * Defines the interface for dynamically loaded annotation plugins.
 * Plugins can provide custom annotation sources.
 */

#ifndef VEP_PLUGIN_HPP
#define VEP_PLUGIN_HPP

#include "annotation_source.hpp"
#include <string>
#include <vector>
#include <memory>
#include <map>

namespace vep {

/**
 * Plugin metadata
 */
struct PluginInfo {
    std::string name;           // Plugin name
    std::string version;        // Version string
    std::string description;    // Description
    std::string author;         // Author/maintainer
    std::string license;        // License
};

/**
 * Plugin configuration
 */
struct PluginConfig {
    std::string plugin_path;    // Path to shared library
    std::string config_string;  // Configuration string passed to plugin
    std::map<std::string, std::string> options;  // Parsed options
};

/**
 * Abstract plugin interface
 *
 * All plugins must implement this interface and export a factory function:
 *   extern "C" VEPPlugin* create_plugin();
 *   extern "C" void destroy_plugin(VEPPlugin* plugin);
 */
class VEPPlugin {
public:
    virtual ~VEPPlugin() = default;

    /**
     * Get plugin information
     */
    virtual PluginInfo get_info() const = 0;

    /**
     * Initialize plugin with configuration
     * @param config Configuration options
     * @return true if initialization successful
     */
    virtual bool initialize(const PluginConfig& config) = 0;

    /**
     * Create annotation sources provided by this plugin
     * @return Vector of annotation sources
     */
    virtual std::vector<std::shared_ptr<AnnotationSource>> create_sources() = 0;

    /**
     * Get required data files
     * @return Map of required file type -> description
     */
    virtual std::map<std::string, std::string> get_required_files() const {
        return {};
    }

    /**
     * Validate configuration
     * @param config Configuration to validate
     * @return Empty string if valid, error message otherwise
     */
    virtual std::string validate_config(const PluginConfig& config) const {
        (void)config;
        return "";
    }
};

/**
 * Plugin factory function types
 */
typedef VEPPlugin* (*CreatePluginFunc)();
typedef void (*DestroyPluginFunc)(VEPPlugin*);

/**
 * Plugin loader - handles dynamic loading of plugins
 */
class PluginLoader {
public:
    PluginLoader() = default;
    ~PluginLoader();

    // Prevent copying
    PluginLoader(const PluginLoader&) = delete;
    PluginLoader& operator=(const PluginLoader&) = delete;

    /**
     * Load a plugin from a shared library
     * @param path Path to .so/.dylib file
     * @param config Configuration string
     * @return true if loaded successfully
     */
    bool load_plugin(const std::string& path, const std::string& config = "");

    /**
     * Add a plugin search directory
     */
    void add_plugin_dir(const std::string& dir);

    /**
     * Load all plugins from registered directories
     */
    void load_plugins_from_dirs();

    /**
     * Get all loaded plugins
     */
    std::vector<VEPPlugin*> get_plugins() const;

    /**
     * Get all annotation sources from all plugins
     */
    std::vector<std::shared_ptr<AnnotationSource>> get_all_sources();

    /**
     * Get plugin count
     */
    size_t plugin_count() const { return plugins_.size(); }

    /**
     * Get error message from last operation
     */
    std::string get_last_error() const { return last_error_; }

private:
    struct LoadedPlugin {
        void* handle;
        VEPPlugin* plugin;
        DestroyPluginFunc destroy_func;
        std::string path;
    };

    std::vector<LoadedPlugin> plugins_;
    std::vector<std::string> plugin_dirs_;
    std::string last_error_;

    void unload_plugin(LoadedPlugin& p);
};

/**
 * Helper macro for plugin implementation
 */
#define VEP_PLUGIN_EXPORT(PluginClass) \
    extern "C" { \
        vep::VEPPlugin* create_plugin() { return new PluginClass(); } \
        void destroy_plugin(vep::VEPPlugin* plugin) { delete plugin; } \
    }

} // namespace vep

#endif // VEP_PLUGIN_HPP
