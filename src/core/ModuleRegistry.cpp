#include "ModuleRegistry.h"
#include "Module.h"

namespace aplaceholder {

ModuleRegistry& ModuleRegistry::instance()
{
    static ModuleRegistry reg;
    return reg;
}

void ModuleRegistry::registerModule(const QString& name, const QString& category, Factory factory)
{
    m_modules[name] = {category, std::move(factory)};
}

std::unique_ptr<Module> ModuleRegistry::create(const QString& name) const
{
    auto it = m_modules.find(name);
    if (it == m_modules.end())
        return nullptr;
    return it->factory();
}

QStringList ModuleRegistry::moduleNames() const
{
    return m_modules.keys();
}

QStringList ModuleRegistry::categories() const
{
    QStringList cats;
    for (const auto& entry : m_modules) {
        if (!cats.contains(entry.category))
            cats.append(entry.category);
    }
    cats.sort();
    return cats;
}

QStringList ModuleRegistry::modulesInCategory(const QString& category) const
{
    QStringList names;
    for (auto it = m_modules.begin(); it != m_modules.end(); ++it) {
        if (it->category == category)
            names.append(it.key());
    }
    names.sort();
    return names;
}

} // namespace aplaceholder
