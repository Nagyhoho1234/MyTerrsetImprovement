#pragma once

#include <QString>
#include <QMap>
#include <memory>
#include <functional>
#include <vector>

namespace aplaceholder {

class Module;

// Central registry of all available modules
class ModuleRegistry {
public:
    using Factory = std::function<std::unique_ptr<Module>()>;

    static ModuleRegistry& instance();

    // Register a module factory
    void registerModule(const QString& name, const QString& category, Factory factory);

    // Create module by name
    std::unique_ptr<Module> create(const QString& name) const;

    // Query
    QStringList moduleNames() const;
    QStringList categories() const;
    QStringList modulesInCategory(const QString& category) const;

    // Auto-registration helper
    template<typename T>
    struct Registrar {
        Registrar() {
            T tmp;
            ModuleRegistry::instance().registerModule(
                tmp.name(), tmp.category(),
                []() { return std::make_unique<T>(); }
            );
        }
    };

private:
    ModuleRegistry() = default;

    struct Entry {
        QString category;
        Factory factory;
    };
    QMap<QString, Entry> m_modules;
};

// Macro for auto-registering modules
#define REGISTER_MODULE(ClassName) \
    static aplaceholder::ModuleRegistry::Registrar<ClassName> s_##ClassName##_registrar;

} // namespace aplaceholder
