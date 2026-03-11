#include "Module.h"

namespace aplaceholder {

void Module::setParameter(const QString& key, const QVariant& value)
{
    m_params[key] = value;
}

QVariant Module::parameter(const QString& key) const
{
    return m_params.value(key);
}

void Module::reportProgress(double progress, const QString& msg)
{
    if (m_progressCb)
        m_progressCb(progress, msg);
}

} // namespace aplaceholder
