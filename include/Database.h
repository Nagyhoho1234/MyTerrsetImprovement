#pragma once

#include <QString>
#include <QStringList>
#include <QVariant>
#include <QMap>
#include <vector>
#include <memory>

struct sqlite3;

namespace aplaceholder {

// Row = map of column_name -> value
using DbRow = QMap<QString, QVariant>;

// Embedded database manager (replaces TerrSet's Access/Jet Engine)
// Zero external dependencies — SQLite is compiled into the application
class Database {
public:
    Database();
    ~Database();

    // Open/create a database file (.db)
    bool open(const QString& path);
    void close();
    bool isOpen() const { return m_db != nullptr; }
    QString path() const { return m_path; }

    // Schema operations
    bool createTable(const QString& name, const QStringList& columns,
                     const QStringList& types);
    bool dropTable(const QString& name);
    QStringList tables() const;
    QStringList columns(const QString& table) const;
    QStringList columnTypes(const QString& table) const;
    int rowCount(const QString& table) const;

    // Data operations
    bool insert(const QString& table, const DbRow& row);
    bool insertBatch(const QString& table, const std::vector<DbRow>& rows);
    bool update(const QString& table, const DbRow& values,
                const QString& whereClause);
    bool deleteRows(const QString& table, const QString& whereClause);

    // Query
    std::vector<DbRow> query(const QString& sql) const;
    std::vector<DbRow> selectAll(const QString& table,
                                  const QString& orderBy = {}) const;
    std::vector<DbRow> filter(const QString& table,
                               const QString& whereClause,
                               const QString& orderBy = {}) const;

    // Raw SQL execution (for advanced users)
    bool execute(const QString& sql);

    // Import/Export
    bool importCSV(const QString& table, const QString& csvPath,
                   bool hasHeader = true, QChar delimiter = ',');
    bool exportCSV(const QString& table, const QString& csvPath,
                   QChar delimiter = ',') const;

    // Link to raster — associate attribute data with raster feature IDs
    bool linkToRaster(const QString& table, const QString& rasterPath,
                      const QString& idField);

    // Calculate field values (like TerrSet's Calculate Field Values)
    bool calculateField(const QString& table, const QString& field,
                        const QString& expression);

    // Error handling
    QString lastError() const { return m_lastError; }

private:
    sqlite3* m_db = nullptr;
    QString m_path;
    mutable QString m_lastError;
};

// Project database — manages the project's central database
// Each project has one .db file containing all attribute tables
class ProjectDatabase {
public:
    static ProjectDatabase& instance();

    bool openProject(const QString& projectDir);
    void closeProject();

    Database& database() { return m_db; }
    const Database& database() const { return m_db; }

    QString projectDir() const { return m_projectDir; }

private:
    ProjectDatabase() = default;
    Database m_db;
    QString m_projectDir;
};

} // namespace aplaceholder
