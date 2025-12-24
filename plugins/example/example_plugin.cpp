/**
 * Example VEP Plugin
 *
 * This is a template for creating custom annotation plugins.
 * Compile with:
 *   g++ -std=c++17 -shared -fPIC -I../../include example_plugin.cpp -o example_plugin.so
 */

#include "plugin.hpp"
#include <sstream>

namespace vep {

/**
 * Example annotation source - adds a custom field to all variants
 */
class ExampleSource : public AnnotationSource {
public:
    ExampleSource(const std::string& message)
        : message_(message) {}

    std::string name() const override { return "example"; }
    std::string type() const override { return "custom"; }
    std::string description() const override {
        return "Example plugin annotation source";
    }

    bool is_ready() const override { return true; }

    void initialize() override {
        // No initialization needed for this example
    }

    void annotate(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt,
        const Transcript* transcript,
        std::map<std::string, std::string>& annotations
    ) override {
        // Add custom annotation
        annotations["example:message"] = message_;
        annotations["example:variant_id"] = chrom + ":" + std::to_string(pos) +
                                           ":" + ref + ":" + alt;

        if (transcript) {
            annotations["example:transcript"] = transcript->id;
        }
    }

    std::vector<std::string> get_fields() const override {
        return {
            "example:message",
            "example:variant_id",
            "example:transcript"
        };
    }

    bool is_thread_safe() const override { return true; }

private:
    std::string message_;
};

/**
 * Example plugin implementation
 */
class ExamplePlugin : public VEPPlugin {
public:
    PluginInfo get_info() const override {
        PluginInfo info;
        info.name = "example";
        info.version = "1.0.0";
        info.description = "Example plugin demonstrating VEP plugin interface";
        info.author = "VEP Development Team";
        info.license = "MIT";
        return info;
    }

    bool initialize(const PluginConfig& config) override {
        // Get message from config, default to "Hello from example plugin"
        auto it = config.options.find("message");
        if (it != config.options.end()) {
            message_ = it->second;
        } else {
            message_ = "Hello from example plugin";
        }
        return true;
    }

    std::vector<std::shared_ptr<AnnotationSource>> create_sources() override {
        return {
            std::make_shared<ExampleSource>(message_)
        };
    }

    std::map<std::string, std::string> get_required_files() const override {
        return {};  // This plugin doesn't require any data files
    }

    std::string validate_config(const PluginConfig& config) const override {
        (void)config;
        return "";  // Any config is valid for this example
    }

private:
    std::string message_;
};

// Export the plugin factory functions
VEP_PLUGIN_EXPORT(ExamplePlugin)

} // namespace vep
