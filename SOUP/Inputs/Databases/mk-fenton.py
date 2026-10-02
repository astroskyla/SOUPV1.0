# This generates the corrected Fenton reaction database, following the Gallard & De Laat (2000) approach.
# Kinetic reactions use the aggregated [Fe(III)] (all the iron(III) forms lumped together), while equilibrium reactions use the individual iron(III) species separately.

import sqlite3
import os
import pandas as pd

COLUMNS = [
    "id",
    "reactants",
    "products",
    "forward_rate_constant",
    "forward_units",
    "forward_reference",
    "reverse_rate_constant",
    "reverse_units",
    "reverse_reference",
    "equilibrium_constant",
    "equilibrium_units",
    "equilibrium_reference",
]

def create_fenton_database(db_name="fenton_corrected"):
    # Purpose: writes the corrected Fenton reaction database, with equilibrium reactions using individual iron(III) species and kinetic reactions using their aggregated total.
    # Inputs are: db_name, the database filename without its .db extension (default 'fenton_corrected').
    # Outputs are: db_path, the full path to the database file created. A SQLite database encoding the corrected Fenton reaction network is written as a side effect, following the Gallard & De Laat approach: the individual iron(III) species (Fe3, FeIIIOH, etc.) are used in the equilibrium reactions, while the aggregated species FeIII_agg (= Fe3 + FeIIIOH + FeIIIOHOH + 2*FeFeIIIOHOH, calculated dynamically during the simulation) is used in the kinetic rate expressions.
    # Called by: the script's entry point, at the bottom of this file (if __name__ == "__main__").

    reactions_raw = [
        # EQUILIBRIUM REACTIONS - use the individual species

        # Iron hydrolysis equilibria (Table 1, reactions I, II, III)
        {
            "reactants": "(1)Fe3 + (1)H2O",
            "products": "(1)FeIIIOH + (1)HPLUS1",
            "equilibrium_constant": 2.9e-3,
            "equilibrium_units": "M",
            "equilibrium_reference": "Milburn and Vosburgh (1955)",
        },
        {
            "reactants": "(1)Fe3 + (2)H2O",
            "products": "(1)FeIIIOHOH + (2)HPLUS1",
            "equilibrium_constant": 7.62e-7,
            "equilibrium_units": "M^2",
            "equilibrium_reference": "Milburn and Vosburgh (1955)",
        },
        {
            "reactants": "(2)Fe3 + (2)H2O",
            "products": "(1)FeFeIIIOHOH + (2)HPLUS1",
            "equilibrium_constant": 0.8e-3,
            "equilibrium_units": "M",
            "equilibrium_reference": "Knight and Sylva (1975)",
        },

        # Fe(III)-peroxo complex formation (Table 1, reactions IV, V)
        {
            "reactants": "(1)Fe3 + (1)H2O2",
            "products": "(1)I1 + (1)HPLUS1",
            "equilibrium_constant": 3.1e-3,
            "equilibrium_units": "dimensionless",
            "equilibrium_reference": "Gallard et al. (1999)",
        },
        {
            "reactants": "(1)FeIIIOH + (1)H2O2",
            "products": "(1)I2 + (1)HPLUS1",
            "equilibrium_constant": 2.0e-4,
            "equilibrium_units": "dimensionless",
            "equilibrium_reference": "Gallard et al. (1999)",
        },

        # KINETIC REACTIONS - use the aggregated species where appropriate

        # Fe(III)-peroxo complex decomposition - rate limiting (Table 1, reaction c)
        {
            "reactants": "(1)I1",
            "products": "(1)FeII + (1)HOOdot",
            "forward_rate_constant": 2.7e-3,
            "forward_units": "s^-1",
            "forward_reference": "De Laat and Gallard (1999)",
        },
        {
            "reactants": "(1)I2",
            "products": "(1)FeII + (1)HOOdot + (1)OHdot",
            "forward_rate_constant": 2.7e-3,
            "forward_units": "s^-1",
            "forward_reference": "De Laat and Gallard (1999)",
        },

        # Classical Fenton reaction (Table 1, reaction o) - uses the aggregated Fe(III)
        {
            "reactants": "(1)FeII + (1)H2O2",
            "products": "(1)FeIII_agg + (1)OHdot + (1)OHMINUS1",
            "forward_rate_constant": 63.0,
            "forward_units": "M^-1s^-1",
            "forward_reference": "Gallard et al. (1998)",
        },

        # Iron-radical reactions (Table 1, reactions 1, 3, 3', 4, 4') - use the aggregated Fe(III)
        {
            "reactants": "(1)FeII + (1)OHdot",
            "products": "(1)FeIII_agg + (1)OHMINUS1",
            "forward_rate_constant": 3.0e8,
            "forward_units": "M^-1s^-1",
            "forward_reference": "Stuglik and Zagorski (1981)",
        },
        {
            "reactants": "(1)FeII + (1)HOOdot",
            "products": "(1)FeIII_agg + (1)H2O2 + (1)OHMINUS1",
            "forward_rate_constant": 1.2e6,
            "forward_units": "M^-1s^-1",
            "forward_reference": "Jayson et al. (1969)",
        },
        {
            "reactants": "(1)FeII + (1)OOdot + (1)HPLUS1",
            "products": "(1)FeIII_agg + (1)H2O2 + (1)OHMINUS1",
            "forward_rate_constant": 1.0e7,
            "forward_units": "M^-1s^-1",
            "forward_reference": "Rush and Bielski (1985)",
        },
        {
            "reactants": "(1)FeIII_agg + (1)HOOdot",
            "products": "(1)FeII + (1)OO + (1)HPLUS1",
            "forward_rate_constant": 1.0e3,
            "forward_units": "M^-1s^-1",
            "forward_reference": "Rush and Bielski (1985)",
        },
        {
            "reactants": "(1)FeIII_agg + (1)OOdot",
            "products": "(1)FeII + (1)OO",
            "forward_rate_constant": 5.0e7,
            "forward_units": "M^-1s^-1",
            "forward_reference": "Rothschild and Allen (1958)",
        },

        # RADICAL REACTIONS - no iron species involved

        {
            "reactants": "(1)OHdot + (1)H2O2",
            "products": "(1)HOOdot + (1)H2O",
            "forward_rate_constant": 3.3e7,
            "forward_units": "M^-1s^-1",
            "forward_reference": "Christensen et al. (1982)",
        },
        {
            "reactants": "(1)HOOdot",
            "products": "(1)OOdot + (1)HPLUS1",
            "forward_rate_constant": 1.58e5,
            "forward_units": "s^-1",
            "forward_reference": "Bielski et al. (1985)",
        },
        {
            "reactants": "(1)OOdot + (1)HPLUS1",
            "products": "(1)HOOdot",
            "forward_rate_constant": 1.0e10,
            "forward_units": "M^-1s^-1",
            "forward_reference": "Bielski et al. (1985)",
        },
        {
            "reactants": "(2)HOOdot",
            "products": "(1)H2O2 + (1)OO",
            "forward_rate_constant": 8.3e5,
            "forward_units": "M^-1s^-1",
            "forward_reference": "Bielski et al. (1985)",
        },
        {
            "reactants": "(1)HOOdot + (1)OOdot + (1)H2O",
            "products": "(1)H2O2 + (1)OO + (1)OHMINUS1",
            "forward_rate_constant": 9.7e7,
            "forward_units": "M^-1s^-1",
            "forward_reference": "Bielski et al. (1985)",
        },
        {
            "reactants": "(1)OHdot + (1)HOOdot",
            "products": "(1)H2O + (1)OO",
            "forward_rate_constant": 7.1e9,
            "forward_units": "M^-1s^-1",
            "forward_reference": "Sehested et al. (1968)",
        },
        {
            "reactants": "(1)OHdot + (1)OOdot",
            "products": "(1)OHMINUS1 + (1)OO",
            "forward_rate_constant": 1.01e10,
            "forward_units": "M^-1s^-1",
            "forward_reference": "Sehested et al. (1968)",
        },
        {
            "reactants": "(2)OHdot",
            "products": "(1)H2O2",
            "forward_rate_constant": 5.2e9,
            "forward_units": "M^-1s^-1",
            "forward_reference": "Sehested et al. (1968)",
        },

        # ATRAZINE OXIDATION REACTIONS

        {
            "reactants": "(1)atrazine + (1)OHdot",
            "products": "(1)Pi",
            "forward_rate_constant": 3.0e9,
            "forward_units": "M^-1s^-1",
            "forward_reference": "Gallard et al. (1998)",
        },
        {
            "reactants": "(1)Pi + (1)OHdot",
            "products": "(1)Pix",
            "forward_rate_constant": 3.0e9,
            "forward_units": "M^-1s^-1",
            "forward_reference": "This study",
        },
        {
            "reactants": "(1)Pix + (1)OHdot",
            "products": "(1)products",
            "forward_rate_constant": 1.0e5,
            "forward_units": "M^-1s^-1",
            "forward_reference": "This study",
        },
    ]

    # Every reaction is normalised to the full column schema here, so any field not given above is filled with None, and each reaction gets a unique id.
    reactions = []
    for i, r in enumerate(reactions_raw, start=1):
        row = {col: None for col in COLUMNS}  # fill all fields with None
        row.update(r)                         # overwrite with provided values
        row["id"] = i                         # assign unique ID
        reactions.append(row)

    df = pd.DataFrame(reactions, columns=COLUMNS)

    db_dir = "../Inputs/Databases"
    os.makedirs(db_dir, exist_ok=True)

    db_path = os.path.join(db_dir, f"{db_name}.db")
    if os.path.exists(db_path):
        os.remove(db_path)

    conn = sqlite3.connect(db_path)
    df.to_sql("reactions", conn, index=False)
    conn.close()

    print(f"Fenton database created: {db_path}")
    print(f"   - {len(df)} reactions")
    print(f"   - Individual species: Fe3, FeIIIOH, FeIIIOHOH, FeFeIIIOHOH, I1, I2")
    print(f"   - Aggregated species: FeIII_agg (calculated dynamically)")
    print(f"   - Equilibrium reactions use individual species")
    print(f"   - Kinetic reactions use aggregated species")
    print(f"   - Follows Gallard & De Laat (2000) approach exactly")

    return db_path


if __name__ == "__main__":
    create_fenton_database()
