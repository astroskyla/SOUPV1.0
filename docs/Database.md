# Database

The SOUP Database is designed to host information relating to reaction networks. To view the master database, create your own custom database, or add reactions to the master database [click here](http://localhost:5000/)

!!! tip "database links"
    The database links on this page point to the Flask database app included in the SOUP repository. They only work if you have cloned the repository and started it locally. See [Running the Model](RunMe.md) for setup instructions. They will not work directly from this hosted documentation site. We apologise for the inconvinience and hope to fix this in the future.

---

## Description of the Database

The SOUP Database is composed of a series of chemical reactions along with their kinetic data, including foward rate constants, reverse rate constants, and equilibrium constants. Where possible, orders of reaction with respect to each reactant are provided, where not possible reaction stoichiometry is included. Each reaction has an associated numerical reference. References can be found here. Alongside this users can find the species name used in the database along with the written names and the associated SMILES by looking here. 

The database provided here is what we term the "master database", this is a culmination of many smaller reaction networks. Examples of which include:

- Iron Hydrolysis (The Fenton Reaction)
- Equilibrium (Comparison to GWB)
- Zymonic Acid Speciation 
- Prebiotic Chemical Network

Where each of these smaller networks can be forked off the master database. Forks can be made by following the procedure outlined below.

---

## Forking the Database 

**Step 1.** Go to view the master database by [clicking here](http://localhost:5000/).

**Step 2.** Next to each reaction id there is a tickbox. Select the reactions you would like to be present in your forked version of the database.

**Step 3.** Once you have selected all the desired reactions scroll to the bottom of the page and select:

    Export Selection

The database file will be downloaded to your computer's default download folder.

**Step 4.** After downloading, please move the file into the repository folder: SOUP/Inputs/Databases

---


## How to contribute

**Step 1.** Go to view the master database by [clicking here](http://localhost:5000/).

**Step 2.** Select the "add new reaction" button above the database at the top of the page.

**Step 3.** Fill in the questions as they appear. The "Units" fields are free text i.e. nothing currently validates them, so it is on you to keep things consistent. We recommend molar units (e.g., M, mM, µM, nM) for concentrations, time units (e.g., s, min, hr) for rate constants, or compound forms like M⁻¹s⁻¹ for bimolecular rates. This matches what is already used throughout the database. Please avoid mass-based units (e.g., mg/mL, % w/v), since nothing downstream converts them.

**Step 4.** Scroll down to the bottom of the questionaire. Here you will have two options. 

1. Save Reaction Locally
2. Submit Reaction for Approval

*Button 1* will save the reaction to your local version of the master database.

*Button 2* will download the added reaction as a .db file and open up your email browser. All you will need to do is attach the downloaded file to this email and click send. From here your reaction will be checked for it's validity and if successful it will be added to the public database. Additions will be made on a rolling basis. 

---