#include "Module.h"
#include "ModuleRegistry.h"

namespace aplaceholder {

class MediaViewerModule : public Module {
public:
    QString name() const override { return "MEDIAVIEWER"; }
    QString description() const override {
        return "Media and animation viewer for raster time series and animations.";
    }
    QString category() const override { return "File"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {};
    }

    bool execute() override {
        reportProgress(1.0, "Media Viewer ready.");
        return true;
    }
};

REGISTER_MODULE(MediaViewerModule)

} // namespace aplaceholder
