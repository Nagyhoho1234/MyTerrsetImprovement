#include <QApplication>
#include <QTabWidget>
#include <iostream>
#include "modelers/LcmWizard.h"
#include "modelers/HbmWizard.h"
#include "modelers/GeOsirisWizard.h"
#include "modelers/EtmWizard.h"

using namespace aplaceholder;

struct WizardTest {
    const char* name;
    int expectedTabs;
    bool passed = false;
};

int main(int argc, char* argv[])
{
    QApplication app(argc, argv);

    int passed = 0, total = 0;

    auto runTest = [&](const char* name, auto* wizard, int expectedTabs) {
        total++;
        auto* tabs = wizard->findChild<QTabWidget*>();
        int tabCount = tabs ? tabs->count() : 0;

        // Collect tab names
        std::string tabNames;
        if (tabs) {
            for (int i = 0; i < tabs->count(); i++) {
                if (i > 0) tabNames += ", ";
                tabNames += tabs->tabText(i).toStdString();
            }
        }

        bool ok = (tabCount == expectedTabs);
        std::cout << total << ". " << name << "... "
                  << (ok ? "PASS" : "FAIL")
                  << " (" << tabCount << " tabs";
        if (expectedTabs != tabCount)
            std::cout << ", expected " << expectedTabs;
        std::cout << ": " << tabNames << ")" << std::endl;

        if (ok) passed++;

        // Test that wizard has a title
        bool hasTitle = !wizard->windowTitle().isEmpty();
        total++;
        std::cout << "   Title: \"" << wizard->windowTitle().toStdString() << "\"... "
                  << (hasTitle ? "PASS" : "FAIL") << std::endl;
        if (hasTitle) passed++;

        // Test that wizard can be shown (doesn't crash)
        total++;
        wizard->show();
        wizard->hide();
        std::cout << "   Show/Hide... PASS" << std::endl;
        passed++;

        delete wizard;
    };

    std::cout << "=== WIZARD TESTS ===" << std::endl;

    runTest("LCM Wizard",      new LcmWizard(),      6);
    runTest("HBM Wizard",      new HbmWizard(),      4);
    runTest("GeOSIRIS Wizard", new GeOsirisWizard(), 7);
    runTest("ETM Wizard",      new EtmWizard(),      3);

    std::cout << "\n=== WIZARD TEST RESULTS ===" << std::endl;
    std::cout << "Passed: " << passed << "/" << total << std::endl;
    std::cout << (passed == total ? "ALL TESTS PASSED" : "SOME TESTS FAILED") << std::endl;

    return (passed == total) ? 0 : 1;
}
