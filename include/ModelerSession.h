#pragma once

#include <QString>
#include <QVariant>
#include <QMap>
#include <memory>

namespace aplaceholder {

class Raster;

/// Shared state container for a modeler wizard session.
/// Tabs communicate through key-value pairs stored here.
class ModelerSession {
public:
    void set(const QString& key, const QVariant& value) { m_params[key] = value; }
    QVariant get(const QString& key, const QVariant& defaultValue = {}) const {
        return m_params.value(key, defaultValue);
    }
    bool contains(const QString& key) const { return m_params.contains(key); }
    void remove(const QString& key) { m_params.remove(key); }

    void setRaster(const QString& key, std::shared_ptr<Raster> raster) {
        m_rasters[key] = std::move(raster);
    }
    std::shared_ptr<Raster> raster(const QString& key) const {
        return m_rasters.value(key);
    }

    void clear() { m_params.clear(); m_rasters.clear(); }

private:
    QMap<QString, QVariant> m_params;
    QMap<QString, std::shared_ptr<Raster>> m_rasters;
};

} // namespace aplaceholder
