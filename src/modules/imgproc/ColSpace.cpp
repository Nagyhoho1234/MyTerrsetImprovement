#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>

namespace aplaceholder {

class ColSpaceModule : public Module {
public:
    QString name() const override { return "COLSPACE"; }
    QString description() const override {
        return "Color space transformation. Converts between RGB and HSI (Hue/Saturation/Intensity) "
               "or HSV (Hue/Saturation/Value) color models.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_bands", "Input bands (comma-separated: R,G,B or H,S,I/V)",
                "Three comma-separated band image file paths"),
            ParameterDef::file("output_prefix", "Output prefix",
                "Prefix for the 3 output band files (e.g., prefix_H, prefix_S, prefix_I)"),
            ParameterDef::combo("direction", "Transformation direction",
                {"RGB_to_HSI", "HSI_to_RGB", "RGB_to_HSV", "HSV_to_RGB"}, 0,
                "Color space conversion direction"),
        };
    }

    bool execute() override {
        // ------------------------------------------------------------------
        // 1. Read input bands
        // ------------------------------------------------------------------
        QString bandsParam = parameter("input_bands").toString();
        QStringList bandPaths = bandsParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : bandPaths) p = p.trimmed();

        if (bandPaths.size() != 3) {
            setError("Exactly 3 input bands are required");
            return false;
        }

        std::vector<std::unique_ptr<Raster>> bandRasters;
        for (const auto& path : bandPaths) {
            auto r = GdalIO::read(path);
            if (!r) {
                setError("Failed to read band image: " + path);
                return false;
            }
            bandRasters.push_back(std::move(r));
        }

        int cols = bandRasters[0]->cols();
        int rows = bandRasters[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = bandRasters[0]->hasNoData();
        double noData = bandRasters[0]->noDataValue();

        for (int b = 1; b < 3; ++b) {
            if (bandRasters[b]->cols() != cols || bandRasters[b]->rows() != rows) {
                setError("All input bands must have the same dimensions");
                return false;
            }
        }

        int direction = parameter("direction").toInt();
        QString prefix = parameter("output_prefix").toString();

        reportProgress(0.1, "Transforming color space...");

        // ------------------------------------------------------------------
        // 2. Create 3 output rasters
        // ------------------------------------------------------------------
        Raster out1(cols, rows, 1, DataType::Float64);
        Raster out2(cols, rows, 1, DataType::Float64);
        Raster out3(cols, rows, 1, DataType::Float64);

        out1.setGeoTransform(bandRasters[0]->geoTransform());
        out1.setProjection(bandRasters[0]->projection());
        out1.setNoDataValue(noData);
        out2.setGeoTransform(bandRasters[0]->geoTransform());
        out2.setProjection(bandRasters[0]->projection());
        out2.setNoDataValue(noData);
        out3.setGeoTransform(bandRasters[0]->geoTransform());
        out3.setProjection(bandRasters[0]->projection());
        out3.setNoDataValue(noData);

        const auto& band1 = bandRasters[0]->data(0);
        const auto& band2 = bandRasters[1]->data(0);
        const auto& band3 = bandRasters[2]->data(0);

        auto& d1 = out1.data(0);
        auto& d2 = out2.data(0);
        auto& d3 = out3.data(0);

        static constexpr double PI = 3.14159265358979323846;
        static constexpr double TWO_PI = 2.0 * PI;

        // ------------------------------------------------------------------
        // 3. Pixel-by-pixel transformation
        // ------------------------------------------------------------------
        QString suffix1, suffix2, suffix3;

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && band1[i] == noData) {
                d1[i] = noData;
                d2[i] = noData;
                d3[i] = noData;
                continue;
            }

            double a = band1[i];
            double b = band2[i];
            double c = band3[i];

            switch (direction) {
            case 0: {
                // RGB to HSI
                double R = a, G = b, B = c;
                double I = (R + G + B) / 3.0;
                double minRGB = std::min({R, G, B});
                double S = (I > 1e-10) ? 1.0 - minRGB / I : 0.0;

                double num = 0.5 * ((R - G) + (R - B));
                double den = std::sqrt((R - G) * (R - G) + (R - B) * (G - B));
                double H;
                if (den < 1e-10) {
                    H = 0.0;
                } else {
                    double theta = std::acos(std::clamp(num / den, -1.0, 1.0));
                    H = (B <= G) ? theta : (TWO_PI - theta);
                }
                // Normalize H to [0, 360]
                H = H * 180.0 / PI;

                d1[i] = H;
                d2[i] = S;
                d3[i] = I;
                suffix1 = "_H"; suffix2 = "_S"; suffix3 = "_I";
                break;
            }
            case 1: {
                // HSI to RGB
                double H = a, S = b, I = c;
                // Convert H from degrees to radians
                H = H * PI / 180.0;

                double R, G, B;
                if (H < TWO_PI / 3.0) {
                    B = I * (1.0 - S);
                    R = I * (1.0 + S * std::cos(H) / std::cos(PI / 3.0 - H));
                    G = 3.0 * I - (R + B);
                } else if (H < 2.0 * TWO_PI / 3.0) {
                    H -= TWO_PI / 3.0;
                    R = I * (1.0 - S);
                    G = I * (1.0 + S * std::cos(H) / std::cos(PI / 3.0 - H));
                    B = 3.0 * I - (R + G);
                } else {
                    H -= 2.0 * TWO_PI / 3.0;
                    G = I * (1.0 - S);
                    B = I * (1.0 + S * std::cos(H) / std::cos(PI / 3.0 - H));
                    R = 3.0 * I - (G + B);
                }

                d1[i] = R;
                d2[i] = G;
                d3[i] = B;
                suffix1 = "_R"; suffix2 = "_G"; suffix3 = "_B";
                break;
            }
            case 2: {
                // RGB to HSV
                double R = a, G = b, B = c;
                double maxVal = std::max({R, G, B});
                double minVal = std::min({R, G, B});
                double delta = maxVal - minVal;

                double V = maxVal;
                double S = (maxVal > 1e-10) ? delta / maxVal : 0.0;
                double H = 0.0;

                if (delta > 1e-10) {
                    if (maxVal == R)
                        H = 60.0 * std::fmod((G - B) / delta + 6.0, 6.0);
                    else if (maxVal == G)
                        H = 60.0 * ((B - R) / delta + 2.0);
                    else
                        H = 60.0 * ((R - G) / delta + 4.0);
                }

                d1[i] = H;
                d2[i] = S;
                d3[i] = V;
                suffix1 = "_H"; suffix2 = "_S"; suffix3 = "_V";
                break;
            }
            case 3: {
                // HSV to RGB
                double H = a, S = b, V = c;
                double C = V * S;
                double Hp = H / 60.0;
                double X = C * (1.0 - std::abs(std::fmod(Hp, 2.0) - 1.0));
                double m = V - C;

                double R1 = 0, G1 = 0, B1 = 0;
                if (Hp < 1.0)       { R1 = C; G1 = X; B1 = 0; }
                else if (Hp < 2.0)  { R1 = X; G1 = C; B1 = 0; }
                else if (Hp < 3.0)  { R1 = 0; G1 = C; B1 = X; }
                else if (Hp < 4.0)  { R1 = 0; G1 = X; B1 = C; }
                else if (Hp < 5.0)  { R1 = X; G1 = 0; B1 = C; }
                else                { R1 = C; G1 = 0; B1 = X; }

                d1[i] = R1 + m;
                d2[i] = G1 + m;
                d3[i] = B1 + m;
                suffix1 = "_R"; suffix2 = "_G"; suffix3 = "_B";
                break;
            }
            }

            if (i % 1000000 == 0)
                reportProgress(0.1 + 0.7 * static_cast<double>(i) / total);
        }

        // Set suffixes based on direction (use final values from loop)
        if (direction == 0) { suffix1 = "_H"; suffix2 = "_S"; suffix3 = "_I"; }
        else if (direction == 1) { suffix1 = "_R"; suffix2 = "_G"; suffix3 = "_B"; }
        else if (direction == 2) { suffix1 = "_H"; suffix2 = "_S"; suffix3 = "_V"; }
        else { suffix1 = "_R"; suffix2 = "_G"; suffix3 = "_B"; }

        // ------------------------------------------------------------------
        // 4. Write 3 output bands
        // ------------------------------------------------------------------
        reportProgress(0.85, "Writing output bands...");

        if (!GdalIO::write(out1, prefix + suffix1 + ".tif")) {
            setError("Failed to write output band 1");
            return false;
        }
        if (!GdalIO::write(out2, prefix + suffix2 + ".tif")) {
            setError("Failed to write output band 2");
            return false;
        }
        if (!GdalIO::write(out3, prefix + suffix3 + ".tif")) {
            setError("Failed to write output band 3");
            return false;
        }

        reportProgress(1.0, "Done.");
        return true;
    }
};

REGISTER_MODULE(ColSpaceModule)

} // namespace aplaceholder
