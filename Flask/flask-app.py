from flask import Flask, render_template_string, send_file, request, redirect, url_for
import sqlite3
import os
import tempfile
import shutil

app = Flask(__name__)

DB_PATH = "master.db"  # your master DB file name

HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>SOUP – Reaction Database</title>

    <!-- Favicon and Logo -->
    <link rel="shortcut icon" href="/static/img/fav2.png">

    <!-- MkDocs Material Base Styles -->
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/material-design-lite/1.3.0/material.indigo-pink.min.css">
    <script defer src="https://code.getmdl.io/1.3.0/material.min.js"></script>
    <link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Roboto:400,700|Roboto+Mono&display=swap">

    <!-- Load Your Custom MkDocs CSS -->
    <link rel="stylesheet" href="/static/css/extra.css">
    <link rel="stylesheet" href="/static/stylesheets/extra.css">

    <style>
        body {
            font-family: "Roboto", sans-serif;
            margin: 0;
            background: #fafafa;
        }
        header {
            background: var(--md-primary-fg-color, #6200ee);
            color: white;
            padding: 16px;
            display: flex;
            align-items: center;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        header img {
            height: 40px;
            margin-right: 10px;
        }
        header h1 {
            font-size: 1.5rem;
            margin: 0;
        }
        .container {
            padding: 20px;
        }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ddd; padding: 8px; }
        th { background-color: #f2f2f2; }
        #searchInput { width: 50%; padding: 8px; margin-bottom: 10px; }
        button {
            background: var(--md-accent-fg-color, #03dac6);
            border: none;
            padding: 10px 20px;
            color: white;
            border-radius: 4px;
            cursor: pointer;
        }
        button:hover { opacity: 0.9; }
        a { color: var(--md-accent-fg-color, #03dac6); }
        .download-tip {
            margin-top: 8px;
            padding: 10px;
            background-color: #eef4f2;
            border: 1px solid #6c9a8b;
            border-radius: 5px;
            font-family: Roboto, sans-serif;
            color: #333;
            max-width: 500px;
        }
    </style>

    <script>
        function searchTable() {
            let input = document.getElementById("searchInput");
            let filter = input.value.toLowerCase();
            let rows = document.querySelectorAll("table tr");
            for (let i = 1; i < rows.length; i++) {
                let cells = rows[i].getElementsByTagName("td");
                let match = false;
                for (let j = 1; j < cells.length; j++) {
                    if (cells[j].textContent.toLowerCase().includes(filter)) {
                        match = true;
                        break;
                    }
                }
                rows[i].style.display = match ? "" : "none";
            }
        }
        function toggleSelectAll(source) {
            checkboxes = document.getElementsByName('reaction_ids');
            for(let i=0; i < checkboxes.length; i++) {
                checkboxes[i].checked = source.checked;
            }
        }
    </script>
</head>
<body>
    <header>
        <img src="/static/img/fav3.png" alt="SOUP Logo">
        <h1>SOUP – Master Database</h1>
    </header>

    <div class="container">
        <input type="text" id="searchInput" onkeyup="searchTable()" placeholder="🔍 Search reactions...">
        <br>
        <a href="{{ url_for('add_reaction_form') }}">➕ Add New Reaction</a>
        <br><br>
        <a href="/download">📥 Download Full Database</a>
        <div class="download-tip">
            The database file will be downloaded to your computer's default download folder (e.g., <em>Downloads</em> folder).<br>
            <strong>After downloading, please move the file into the repository folder:</strong> <code>SOUP/Inputs/Databases</code>
        </div>

        <form method="POST" action="/export">
        <table>
            <tr>
                <th><input type="checkbox" onClick="toggleSelectAll(this)" title="Select/Deselect All"></th>
                {% for header in headers %}
                    <th>{{ header }}</th>
                {% endfor %}
            </tr>
            {% for row in rows %}
                <tr>
                    <td><input type="checkbox" name="reaction_ids" value="{{ row[0] }}"></td>
                    {% for col in row %}
                        <td>{{ col }}</td>
                    {% endfor %}
                </tr>
            {% endfor %}
        </table>
        <br>
        <button type="submit">Export Selection</button>
        </form>
    </div>
</body>
</html>
"""

ADD_REACTION_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Add Reaction - SOUP</title>
    <link rel="shortcut icon" href="/static/img/fav2.png">
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/material-design-lite/1.3.0/material.indigo-pink.min.css">
    <script defer src="https://code.getmdl.io/1.3.0/material.min.js"></script>
    <link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Roboto:400,700|Roboto+Mono&display=swap">
    <link rel="stylesheet" href="/static/css/extra.css">
    <link rel="stylesheet" href="/static/stylesheets/extra.css">
    <style>
        body { font-family: "Roboto", sans-serif; margin: 0; background: #fafafa; }
        header { background: var(--md-primary-fg-color, #6200ee); color: white; padding: 16px; display: flex; align-items: center; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        header img { height: 40px; margin-right: 10px; }
        header h1 { font-size: 1.5rem; margin: 0; }
        .container { padding: 20px; }
        form { max-width: 650px; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }
        label { display: block; margin-top: 15px; font-weight: bold; }
        input[type=text], input[type=number] { width: 100%; padding: 8px; margin-top: 5px; border: 1px solid #ccc; border-radius: 4px; box-sizing: border-box; }
        small { display: block; color: #666; font-size: 12px; margin-top: 4px; }
        button, .email-btn {
            margin-top: 20px;
            background: var(--md-accent-fg-color, #03dac6);
            border: none;
            padding: 10px 20px;
            color: white;
            border-radius: 4px;
            cursor: pointer;
            font-size: 1rem;
            text-decoration: none;
            display: inline-block;
        }
        button:hover, .email-btn:hover { opacity: 0.9; }
        a { color: var(--md-accent-fg-color, #03dac6); display: inline-block; margin-top: 20px; }
    </style>
</head>
<body>
    <header>
        <img src="/static/img/fav3.png" alt="SOUP Logo">
        <h1>SOUP – Add New Reaction</h1>
    </header>

    <div class="container">
        <h2>Add New Reaction (Local Save)</h2>
        <form method="POST" action="/add_reaction">
            
            <label>Reactants (required):
                <input name="reactants" type="text" required placeholder="e.g. (2)H2 + (1)O2">
            </label>
            <small>If known, include the order of reaction in brackets; otherwise use stoichiometry.</small>

            <label>Products (required):
                <input name="products" type="text" required placeholder="e.g. (2)H2O">
            </label>
            <small>If known, include the order of reaction for each reactant (if reversible); otherwise use stoichiometry.</small>

            <label>Forward Rate Constant:
                <input name="forward_rate_constant" type="number" step="any" placeholder="e.g. 1.23-1">
            </label>

            <label>Forward Units:
                <input name="forward_units" type="text" placeholder="Units must be in the form Molar/time, e.g. M/yr">
            </label>

            <label>Forward Reference:
                <input name="forward_reference" type="text" placeholder="e.g. Smith et al. 2021">
            </label>

            <label>Reverse Rate Constant:
                <input name="reverse_rate_constant" type="number" step="any" placeholder="e.g. 4.56e-2">
            </label>

            <label>Reverse Units:
                <input name="reverse_units" type="text" placeholder="Units must be in the form Molar/time, e.g. M/yr">
            </label>

            <label>Reverse Reference:
                <input name="reverse_reference" type="text" placeholder="e.g. Johnson et al. 2020">
            </label>

            <label>Equilibrium Constant:
                <input name="equilibrium_constant" type="number" step="any" placeholder="7.89e-3">
            </label>

            <label>Equilibrium Units:
                <input name="equilibrium_units" type="text" placeholder="Units must be in the form Molar/time, e.g. M/yr. Note equilibrium constants can be dimensionless">
            </label>

            <label>Equilibrium Reference:
                <input name="equilibrium_reference" type="text" placeholder="e.g. Lee et al. 2019">
            </label>

            <button type="submit" name="action" value="local">💾 Save Reaction Locally</button>
            <button type="button" id="exportAndEmailBtn">📥 Submit Reaction for Approval</button>

            <script>
            document.getElementById("exportAndEmailBtn").addEventListener("click", function() {
                const form = document.querySelector("form");
                form.action = "/add_reaction_and_export";
                form.method = "POST";
                form.target = "_self"; // ensures download works
                form.submit();

                // Wait for download to trigger, then open mail client
                setTimeout(() => {
                    window.location.href = "mailto:soup.database.approvals@gmail.com?subject=SOUP%20Reaction%20Submission&body=Please%20attach%20the%20downloaded%20custom_network.db%20file%20to%20this%20email.";
                }, 1500);
            });
            </script>


        </form>

        <p><a href="{{ url_for('index') }}">⬅ Back to Master Database</a></p>
    </div>
</body>
</html>
"""

@app.route("/")
def index():
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("SELECT * FROM reactions")
    rows = c.fetchall()
    headers = [description[0] for description in c.description]
    conn.close()
    return render_template_string(HTML_TEMPLATE, headers=headers, rows=rows)

@app.route("/download")
def download():
    return send_file(DB_PATH, as_attachment=True)

@app.route("/export", methods=["POST"])
def export():
    selected_ids = request.form.getlist("reaction_ids")
    if not selected_ids:
        return redirect(url_for("index"))

    tmp_dir = tempfile.mkdtemp()
    new_db_path = os.path.join(tmp_dir, "custom_network.db")

    master_conn = sqlite3.connect(DB_PATH)
    master_c = master_conn.cursor()

    new_conn = sqlite3.connect(new_db_path)
    new_c = new_conn.cursor()

    new_c.execute("""
    CREATE TABLE reactions (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        reactants TEXT NOT NULL,
        products TEXT NOT NULL,
        forward_rate_constant REAL,
        forward_units TEXT,
        forward_reference TEXT,
        reverse_rate_constant REAL,
        reverse_units TEXT,
        reverse_reference TEXT,
        equilibrium_constant REAL,
        equilibrium_units TEXT,
        equilibrium_reference TEXT
    );
    """)
    new_conn.commit()

    placeholder = ",".join("?" for _ in selected_ids)
    query = f"SELECT * FROM reactions WHERE id IN ({placeholder})"
    master_c.execute(query, selected_ids)
    rows = master_c.fetchall()

    new_c.executemany("""
    INSERT INTO reactions (
        id, reactants, products, forward_rate_constant, forward_units, forward_reference,
        reverse_rate_constant, reverse_units, reverse_reference,
        equilibrium_constant, equilibrium_units, equilibrium_reference
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, rows)
    new_conn.commit()

    master_conn.close()
    new_conn.close()

    response = send_file(new_db_path, as_attachment=True)
    shutil.rmtree(tmp_dir, ignore_errors=True)
    return response

@app.route("/add_reaction", methods=["GET"])
def add_reaction_form():
    return render_template_string(ADD_REACTION_TEMPLATE)

@app.route("/add_reaction", methods=["POST"])
def add_reaction_submit():
    data = {field: request.form.get(field) for field in [
        "reactants", "products", "forward_rate_constant", "forward_units",
        "forward_reference", "reverse_rate_constant", "reverse_units",
        "reverse_reference", "equilibrium_constant", "equilibrium_units",
        "equilibrium_reference"
    ]}

    def parse_float(value):
        try:
            return float(value) if value else None
        except ValueError:
            return None

    for key in ["forward_rate_constant", "reverse_rate_constant", "equilibrium_constant"]:
        data[key] = parse_float(data[key])

    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("""
        INSERT INTO reactions (
            reactants, products,
            forward_rate_constant, forward_units, forward_reference,
            reverse_rate_constant, reverse_units, reverse_reference,
            equilibrium_constant, equilibrium_units, equilibrium_reference
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, tuple(data.values()))
    conn.commit()
    conn.close()

    return redirect(url_for("index"))

@app.route("/add_reaction_and_export", methods=["POST"])
def add_reaction_and_export():
    data = {field: request.form.get(field) for field in [
        "reactants", "products", "forward_rate_constant", "forward_units",
        "forward_reference", "reverse_rate_constant", "reverse_units",
        "reverse_reference", "equilibrium_constant", "equilibrium_units",
        "equilibrium_reference"
    ]}

    def parse_float(value):
        try:
            return float(value) if value else None
        except ValueError:
            return None

    for key in ["forward_rate_constant", "reverse_rate_constant", "equilibrium_constant"]:
        data[key] = parse_float(data[key])

    # Create temporary database
    tmp_dir = tempfile.mkdtemp()
    custom_db_path = os.path.join(tmp_dir, "custom_network.db")

    conn = sqlite3.connect(custom_db_path)
    c = conn.cursor()

    # Create table
    c.execute("""
    CREATE TABLE reactions (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        reactants TEXT NOT NULL,
        products TEXT NOT NULL,
        forward_rate_constant REAL,
        forward_units TEXT,
        forward_reference TEXT,
        reverse_rate_constant REAL,
        reverse_units TEXT,
        reverse_reference TEXT,
        equilibrium_constant REAL,
        equilibrium_units TEXT,
        equilibrium_reference TEXT
    );
    """)

    # Insert single reaction
    c.execute("""
        INSERT INTO reactions (
            reactants, products, forward_rate_constant, forward_units, forward_reference,
            reverse_rate_constant, reverse_units, reverse_reference,
            equilibrium_constant, equilibrium_units, equilibrium_reference
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, tuple(data.values()))

    conn.commit()
    conn.close()

    # Send file and cleanup
    response = send_file(custom_db_path, as_attachment=True)
    shutil.rmtree(tmp_dir, ignore_errors=True)
    return response


if __name__ == "__main__":
    app.run(debug=False, use_reloader=False)
