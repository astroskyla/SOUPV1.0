# Species List

The table below indicates species that are currently included in the database. The species column shows how species are currently shown in the database, this is how users users should add these species when adding reactions to the master database to ensure no duplicates. We hope to have functionaility to detect duplicates in the future, however, this is currently not supported. The SMILES column can be used to make cross comparisons with PubChem.

---

<div id="sheet-container">Loading spreadsheet...</div>

<!-- Load SheetJS -->
<script src="https://cdn.jsdelivr.net/npm/xlsx@0.18.5/dist/xlsx.full.min.js"></script>

<!-- Load DataTables -->
<link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
<script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
<script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

<script>
fetch("../species.xlsx")
  .then(response => {
    if (!response.ok) throw new Error("HTTP " + response.status);
    return response.arrayBuffer();
  })
  .then(data => {
    const workbook = XLSX.read(data, { type: "array" });
    const sheet = workbook.Sheets[workbook.SheetNames[0]];
    
    const jsonData = XLSX.utils.sheet_to_json(sheet, { defval: "" });
    
    if (jsonData.length === 0) {
      document.getElementById("sheet-container").innerText = "Spreadsheet is empty.";
      return;
    }

    const table = document.createElement("table");
    table.setAttribute("id", "excel-table");
    table.classList.add("display");

    const thead = document.createElement("thead");
    const headerRow = document.createElement("tr");
    Object.keys(jsonData[0]).forEach(key => {
      const th = document.createElement("th");
      th.textContent = key;
      headerRow.appendChild(th);
    });
    thead.appendChild(headerRow);
    table.appendChild(thead);

    const tbody = document.createElement("tbody");
    jsonData.forEach(row => {
      const tr = document.createElement("tr");
      Object.values(row).forEach(cellValue => {
        const td = document.createElement("td");
        td.textContent = cellValue;
        tr.appendChild(td);
      });
      tbody.appendChild(tr);
    });
    table.appendChild(tbody);

    const container = document.getElementById("sheet-container");
    container.innerHTML = "";
    container.appendChild(table);

    $(document).ready(function() {
      $("#excel-table").DataTable({
        searching: true,
        paging: true,
        info: true
      });
    });
  })
  .catch(err => {
    document.getElementById("sheet-container").innerText = "Error loading spreadsheet.";
    console.error("Failed to load Excel:", err);
  });
</script>

<style>
#excel-table {
  width: 100% !important;
  border-collapse: collapse;
  font-family: 'Segoe UI', Tahoma, sans-serif;
  font-size: 14px;
  background: #fff;
}

#excel-table th {
  background: #6c9a8b;
  color: white;
  font-weight: 600;
  padding: 8px;
}

#excel-table td {
  padding: 8px;
  border: 1px solid #ddd;
}

#excel-table tr:nth-child(even) {
  background: #f9f9f9;
}

.dataTables_wrapper .dataTables_filter input {
  border-radius: 5px;
  border: 1px solid #ccc;
  padding: 5px;
  margin-bottom: 15px;
}

}

.dataTables_wrapper .dataTables_paginate .paginate_button {
  border-radius: 3px;
  padding: 4px 8px;
  margin: 2px;
}
</style>


<style>
#excel-table {
  width: 100% !important;
  border-collapse: collapse;
  font-family: 'Segoe UI', Tahoma, sans-serif;
  font-size: 14px;
  background: #fff;
}

#excel-table th {
  background: #6c9a8b;
  color: white;
  font-weight: 600;
  padding: 8px;
}

#excel-table td {
  padding: 8px;
  border: 1px solid #ddd;
}

#excel-table tr:nth-child(even) {
  background: #f9f9f9;
}

.dataTables_wrapper .dataTables_filter input {
  border-radius: 5px;
  border: 1px solid #ccc;
  padding: 5px;
}

.dataTables_wrapper .dataTables_paginate .paginate_button {
  border-radius: 3px;
  padding: 4px 8px;
  margin: 2px;
}
</style>