#include "Database.h"
#include <sqlite3.h>
#include <QFile>
#include <QTextStream>
#include <QFileInfo>
#include <QDir>

namespace aplaceholder {

Database::Database() = default;

Database::~Database()
{
    close();
}

bool Database::open(const QString& path)
{
    close();
    int rc = sqlite3_open(path.toUtf8().constData(), &m_db);
    if (rc != SQLITE_OK) {
        m_lastError = QString("Cannot open database: %1")
            .arg(sqlite3_errmsg(m_db));
        sqlite3_close(m_db);
        m_db = nullptr;
        return false;
    }
    m_path = path;

    // Enable WAL mode for better concurrent performance
    execute("PRAGMA journal_mode=WAL");
    // Enable foreign keys
    execute("PRAGMA foreign_keys=ON");
    return true;
}

void Database::close()
{
    if (m_db) {
        sqlite3_close(m_db);
        m_db = nullptr;
    }
    m_path.clear();
}

bool Database::createTable(const QString& name, const QStringList& columns,
                            const QStringList& types)
{
    if (columns.size() != types.size()) {
        m_lastError = "Column count doesn't match type count";
        return false;
    }

    QStringList defs;
    for (int i = 0; i < columns.size(); ++i) {
        defs.append(QString("\"%1\" %2").arg(columns[i], types[i]));
    }

    QString sql = QString("CREATE TABLE IF NOT EXISTS \"%1\" (%2)")
        .arg(name, defs.join(", "));
    return execute(sql);
}

bool Database::dropTable(const QString& name)
{
    return execute(QString("DROP TABLE IF EXISTS \"%1\"").arg(name));
}

QStringList Database::tables() const
{
    QStringList result;
    auto rows = query("SELECT name FROM sqlite_master WHERE type='table' "
                      "AND name NOT LIKE 'sqlite_%' ORDER BY name");
    for (const auto& row : rows) {
        result.append(row["name"].toString());
    }
    return result;
}

QStringList Database::columns(const QString& table) const
{
    QStringList result;
    auto rows = query(QString("PRAGMA table_info(\"%1\")").arg(table));
    for (const auto& row : rows) {
        result.append(row["name"].toString());
    }
    return result;
}

QStringList Database::columnTypes(const QString& table) const
{
    QStringList result;
    auto rows = query(QString("PRAGMA table_info(\"%1\")").arg(table));
    for (const auto& row : rows) {
        result.append(row["type"].toString());
    }
    return result;
}

int Database::rowCount(const QString& table) const
{
    auto rows = query(QString("SELECT COUNT(*) AS cnt FROM \"%1\"").arg(table));
    if (!rows.empty())
        return rows[0]["cnt"].toInt();
    return 0;
}

bool Database::insert(const QString& table, const DbRow& row)
{
    QStringList cols, placeholders;
    QList<QVariant> values;
    for (auto it = row.begin(); it != row.end(); ++it) {
        cols.append(QString("\"%1\"").arg(it.key()));
        placeholders.append("?");
        values.append(it.value());
    }

    QString sql = QString("INSERT INTO \"%1\" (%2) VALUES (%3)")
        .arg(table, cols.join(","), placeholders.join(","));

    sqlite3_stmt* stmt = nullptr;
    int rc = sqlite3_prepare_v2(m_db, sql.toUtf8().constData(), -1, &stmt, nullptr);
    if (rc != SQLITE_OK) {
        m_lastError = sqlite3_errmsg(m_db);
        return false;
    }

    for (int i = 0; i < values.size(); ++i) {
        const auto& v = values[i];
        if (v.isNull()) {
            sqlite3_bind_null(stmt, i + 1);
        } else if (v.typeId() == QMetaType::Int || v.typeId() == QMetaType::LongLong) {
            sqlite3_bind_int64(stmt, i + 1, v.toLongLong());
        } else if (v.typeId() == QMetaType::Double || v.typeId() == QMetaType::Float) {
            sqlite3_bind_double(stmt, i + 1, v.toDouble());
        } else {
            QByteArray ba = v.toString().toUtf8();
            sqlite3_bind_text(stmt, i + 1, ba.constData(), ba.size(), SQLITE_TRANSIENT);
        }
    }

    rc = sqlite3_step(stmt);
    sqlite3_finalize(stmt);

    if (rc != SQLITE_DONE) {
        m_lastError = sqlite3_errmsg(m_db);
        return false;
    }
    return true;
}

bool Database::insertBatch(const QString& table, const std::vector<DbRow>& rows)
{
    if (rows.empty()) return true;

    execute("BEGIN TRANSACTION");
    for (const auto& row : rows) {
        if (!insert(table, row)) {
            execute("ROLLBACK");
            return false;
        }
    }
    return execute("COMMIT");
}

bool Database::update(const QString& table, const DbRow& values,
                       const QString& whereClause)
{
    QStringList sets;
    for (auto it = values.begin(); it != values.end(); ++it) {
        if (it.value().typeId() == QMetaType::QString)
            sets.append(QString("\"%1\"='%2'").arg(it.key(), it.value().toString()));
        else
            sets.append(QString("\"%1\"=%2").arg(it.key(), it.value().toString()));
    }

    QString sql = QString("UPDATE \"%1\" SET %2").arg(table, sets.join(","));
    if (!whereClause.isEmpty())
        sql += " WHERE " + whereClause;
    return execute(sql);
}

bool Database::deleteRows(const QString& table, const QString& whereClause)
{
    QString sql = QString("DELETE FROM \"%1\"").arg(table);
    if (!whereClause.isEmpty())
        sql += " WHERE " + whereClause;
    return execute(sql);
}

std::vector<DbRow> Database::query(const QString& sql) const
{
    std::vector<DbRow> result;
    sqlite3_stmt* stmt = nullptr;

    int rc = sqlite3_prepare_v2(m_db, sql.toUtf8().constData(), -1, &stmt, nullptr);
    if (rc != SQLITE_OK) {
        m_lastError = sqlite3_errmsg(m_db);
        return result;
    }

    int colCount = sqlite3_column_count(stmt);
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        DbRow row;
        for (int i = 0; i < colCount; ++i) {
            QString colName = QString::fromUtf8(sqlite3_column_name(stmt, i));
            int type = sqlite3_column_type(stmt, i);
            switch (type) {
            case SQLITE_INTEGER:
                row[colName] = static_cast<qlonglong>(sqlite3_column_int64(stmt, i));
                break;
            case SQLITE_FLOAT:
                row[colName] = sqlite3_column_double(stmt, i);
                break;
            case SQLITE_TEXT:
                row[colName] = QString::fromUtf8(
                    reinterpret_cast<const char*>(sqlite3_column_text(stmt, i)));
                break;
            case SQLITE_NULL:
            default:
                row[colName] = QVariant();
                break;
            }
        }
        result.push_back(std::move(row));
    }

    sqlite3_finalize(stmt);
    return result;
}

std::vector<DbRow> Database::selectAll(const QString& table,
                                        const QString& orderBy) const
{
    QString sql = QString("SELECT * FROM \"%1\"").arg(table);
    if (!orderBy.isEmpty())
        sql += " ORDER BY " + orderBy;
    return query(sql);
}

std::vector<DbRow> Database::filter(const QString& table,
                                     const QString& whereClause,
                                     const QString& orderBy) const
{
    QString sql = QString("SELECT * FROM \"%1\" WHERE %2").arg(table, whereClause);
    if (!orderBy.isEmpty())
        sql += " ORDER BY " + orderBy;
    return query(sql);
}

bool Database::execute(const QString& sql)
{
    char* errMsg = nullptr;
    int rc = sqlite3_exec(m_db, sql.toUtf8().constData(), nullptr, nullptr, &errMsg);
    if (rc != SQLITE_OK) {
        m_lastError = errMsg ? QString::fromUtf8(errMsg) : "Unknown error";
        sqlite3_free(errMsg);
        return false;
    }
    return true;
}

bool Database::importCSV(const QString& table, const QString& csvPath,
                          bool hasHeader, QChar delimiter)
{
    QFile file(csvPath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        m_lastError = "Cannot open CSV file: " + csvPath;
        return false;
    }

    QTextStream in(&file);
    QStringList headers;

    // Read header or generate column names
    if (hasHeader) {
        QString headerLine = in.readLine();
        headers = headerLine.split(delimiter);
        for (auto& h : headers)
            h = h.trimmed().remove('"');
    }

    // Read first data line to determine types and column count
    QString firstLine = in.readLine();
    if (firstLine.isEmpty()) {
        m_lastError = "CSV file is empty";
        return false;
    }

    QStringList firstValues = firstLine.split(delimiter);
    int colCount = firstValues.size();

    if (!hasHeader) {
        for (int i = 0; i < colCount; ++i)
            headers.append(QString("COL%1").arg(i + 1));
    }

    // Determine types from first row
    QStringList types;
    for (const auto& val : firstValues) {
        QString v = val.trimmed().remove('"');
        bool isInt = false, isDouble = false;
        v.toLongLong(&isInt);
        if (!isInt) v.toDouble(&isDouble);

        if (isInt) types.append("INTEGER");
        else if (isDouble) types.append("REAL");
        else types.append("TEXT");
    }

    // Create table
    dropTable(table);
    if (!createTable(table, headers, types))
        return false;

    // Insert data
    execute("BEGIN TRANSACTION");

    // Insert first row
    auto parseRow = [&](const QStringList& values) {
        DbRow row;
        for (int i = 0; i < qMin(values.size(), headers.size()); ++i) {
            QString v = values[i].trimmed().remove('"');
            if (types[i] == "INTEGER")
                row[headers[i]] = v.toLongLong();
            else if (types[i] == "REAL")
                row[headers[i]] = v.toDouble();
            else
                row[headers[i]] = v;
        }
        insert(table, row);
    };

    parseRow(firstValues);

    // Insert remaining rows
    while (!in.atEnd()) {
        QString line = in.readLine();
        if (line.trimmed().isEmpty()) continue;
        parseRow(line.split(delimiter));
    }

    execute("COMMIT");
    return true;
}

bool Database::exportCSV(const QString& table, const QString& csvPath,
                          QChar delimiter) const
{
    QFile file(csvPath);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        m_lastError = "Cannot create CSV file: " + csvPath;
        return false;
    }

    QTextStream out(&file);
    auto cols = columns(table);

    // Header
    out << cols.join(delimiter) << "\n";

    // Data
    auto rows = selectAll(table);
    for (const auto& row : rows) {
        QStringList vals;
        for (const auto& col : cols) {
            QVariant v = row.value(col);
            if (v.typeId() == QMetaType::QString)
                vals.append(QString("\"%1\"").arg(v.toString()));
            else
                vals.append(v.toString());
        }
        out << vals.join(delimiter) << "\n";
    }

    return true;
}

bool Database::linkToRaster(const QString& table, const QString& rasterPath,
                             const QString& idField)
{
    // Store the link in a metadata table
    execute("CREATE TABLE IF NOT EXISTS _raster_links ("
            "table_name TEXT, raster_path TEXT, id_field TEXT, "
            "PRIMARY KEY(table_name))");

    execute(QString("DELETE FROM _raster_links WHERE table_name='%1'").arg(table));

    DbRow link;
    link["table_name"] = table;
    link["raster_path"] = rasterPath;
    link["id_field"] = idField;
    return insert("_raster_links", link);
}

bool Database::calculateField(const QString& table, const QString& field,
                               const QString& expression)
{
    // SQLite supports expressions directly in UPDATE
    QString sql = QString("UPDATE \"%1\" SET \"%2\" = %3")
        .arg(table, field, expression);
    return execute(sql);
}

// --- ProjectDatabase ---

ProjectDatabase& ProjectDatabase::instance()
{
    static ProjectDatabase inst;
    return inst;
}

bool ProjectDatabase::openProject(const QString& projectDir)
{
    m_projectDir = projectDir;
    QDir dir(projectDir);
    if (!dir.exists())
        dir.mkpath(".");

    QString dbPath = dir.filePath("project.db");
    return m_db.open(dbPath);
}

void ProjectDatabase::closeProject()
{
    m_db.close();
    m_projectDir.clear();
}

} // namespace aplaceholder
