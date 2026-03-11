#pragma once

#include <QString>
#include <QVariant>
#include <QStringList>

namespace aplaceholder {

// Describes a module parameter for auto-generating GUI dialogs
struct ParameterDef {
    enum Type {
        File,           // File picker
        OutputFile,     // Output file picker
        Integer,        // Spin box
        Double,         // Double spin box
        String,         // Text field
        Combo,          // Dropdown
        Bool,           // Checkbox
        RasterLayer,    // Raster layer from project
        MultiRaster,    // Multiple raster layers
        BandSelect,     // Band selector
        Separator       // Visual separator in dialog
    };

    QString key;            // Internal parameter name
    QString label;          // Display label
    Type type;
    QVariant defaultValue;
    QStringList options;    // For Combo type
    double minValue = 0;    // For numeric types
    double maxValue = 999999;
    QString tooltip;
    bool required = true;

    static ParameterDef file(const QString& key, const QString& label,
                              const QString& tooltip = {}) {
        return {key, label, File, {}, {}, 0, 0, tooltip, true};
    }

    static ParameterDef output(const QString& key, const QString& label,
                                const QString& tooltip = {}) {
        return {key, label, OutputFile, {}, {}, 0, 0, tooltip, true};
    }

    static ParameterDef integer(const QString& key, const QString& label,
                                 int defaultVal, int min = 0, int max = 999999,
                                 const QString& tooltip = {}) {
        return {key, label, Integer, defaultVal, {}, (double)min, (double)max, tooltip, true};
    }

    static ParameterDef real(const QString& key, const QString& label,
                              double defaultVal, double min = 0, double max = 999999,
                              const QString& tooltip = {}) {
        return {key, label, Double, defaultVal, {}, min, max, tooltip, true};
    }

    static ParameterDef combo(const QString& key, const QString& label,
                               const QStringList& options, int defaultIndex = 0,
                               const QString& tooltip = {}) {
        return {key, label, Combo, defaultIndex, options, 0, 0, tooltip, true};
    }

    static ParameterDef boolean(const QString& key, const QString& label,
                                 bool defaultVal = false, const QString& tooltip = {}) {
        return {key, label, Bool, defaultVal, {}, 0, 0, tooltip, false};
    }

    static ParameterDef raster(const QString& key, const QString& label,
                                const QString& tooltip = {}) {
        return {key, label, RasterLayer, {}, {}, 0, 0, tooltip, true};
    }

    static ParameterDef separator(const QString& label = {}) {
        return {{}, label, Separator, {}, {}, 0, 0, {}, false};
    }
};

} // namespace aplaceholder
