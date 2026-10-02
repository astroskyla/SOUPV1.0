# This script creates a readable HTML report of reaction fluxes over time, showing how fast each reaction is running at a handful of points during a simulation. Note Claude code assisted in this HTML report owing to my own lack of ability to write HTML well.

import os
import numpy as np
from datetime import datetime
import re

def format_species_name(name):
    # Purpose: converts a species name from the internal database format to the display format.
    # Inputs are: name, a species name in the internal database format, e.g. 'HPLUS1'.
    # Outputs are: the same name converted to the display format, e.g. 'H+1'.
    # Called by: not called anywhere in the repo (ReactionRatesLogger uses its own identical method instead).
    return name.replace("PLUS", "+").replace("MINUS", "-")

class ReactionRatesLogger:
    # This logs reaction rates at a limited number of points during a simulation, and writes them out as a formatted HTML report once the run finishes.
    # Instantiated by: UniversalRobustODESystem.__init__(), in differentials.py; create_simple_fenton_ode_system(), in main.py.

    def __init__(self, database_species, differential_equations, rate_laws, output_file='./Outputs/reaction_rates.html'):
        # Purpose: sets up a new logger and starts the HTML report file.
        # Inputs are: database_species, the list of species names; differential_equations, the equations built by construct_differential_equations(); rate_laws, the rate-law list, used to identify which reaction each term in an equation comes from; output_file, where the HTML report is written.
        # Outputs are: none directly; the output HTML file is created (any old file at that path is removed first) and given its opening template.
        # Called by: automatically, when a ReactionRatesLogger is created - see "Instantiated by" above.
        self.database_species = database_species
        self.differential_equations = differential_equations
        self.rate_laws = rate_laws
        self.output_file = output_file
        self.species_to_index = {species: i for i, species in enumerate(database_species)}

        self.setup_html_file()

        # This tracks which timepoints have already been logged, so should_log_time() doesn't log the same point twice.
        self.logged_times = set()

    def setup_html_file(self):
        # Purpose: writes the report's opening HTML (page header and CSS/JS styling).
        # Inputs are: none.
        # Outputs are: none directly; the output HTML file is (re)created with its page header and styling written to it.
        # Called by: __init__(), in this file.

        os.makedirs('./Outputs', exist_ok=True)

        if os.path.exists(self.output_file):
            os.remove(self.output_file)

        with open(self.output_file, 'w') as f:
            f.write(f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Chemical Reaction Rates Analysis</title>
    <style>
        :root {{
            --primary-light: #DDE7D8;
            --primary-medium: #6C9A8B;
            --primary-dark: #173F3F;
            --accent-brown: #5C4033;
            --accent-orange: #F37B3A;
            --production-green: #27ae60;
            --loss-red: #e74c3c;
        }}

        body {{
            font-family: 'Segoe UI', -apple-system, BlinkMacSystemFont, sans-serif;
            line-height: 1.6;
            color: var(--primary-dark);
            background: linear-gradient(135deg, var(--primary-light) 0%, #f0f5f1 100%);
            margin: 0;
            padding: 20px;
        }}

        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(92, 64, 51, 0.1);
            padding: 30px;
        }}

        .header {{
            text-align: center;
            margin-bottom: 40px;
            padding-bottom: 20px;
            border-bottom: 3px solid var(--primary-medium);
        }}

        h1 {{
            color: var(--primary-dark);
            font-size: 2.5em;
            margin-bottom: 10px;
        }}

        .subtitle {{
            color: var(--accent-brown);
            font-size: 1.1em;
        }}

        .time-section {{
            margin: 40px 0;
            border: 1px solid var(--primary-light);
            border-radius: 10px;
            overflow: hidden;
            box-shadow: 0 4px 6px rgba(23, 63, 63, 0.1);
        }}

        .time-header {{
            background: linear-gradient(135deg, var(--primary-medium), var(--primary-dark));
            color: white;
            padding: 15px 25px;
            font-size: 1.3em;
            font-weight: 600;
        }}

        .species-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(350px, 1fr));
            gap: 20px;
            padding: 20px;
            background: #f8faf9;
        }}

        .species-card {{
            background: white;
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 2px 8px rgba(23, 63, 63, 0.08);
            border-left: 4px solid var(--primary-medium);
        }}

        .species-name {{
            font-size: 1.2em;
            font-weight: 600;
            color: var(--primary-dark);
            margin-bottom: 10px;
        }}

        .rate-info {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 15px;
        }}

        .concentration {{
            color: var(--accent-brown);
            font-size: 0.9em;
        }}

        .net-rate {{
            font-weight: 600;
            font-size: 1.1em;
            padding: 5px 12px;
            border-radius: 20px;
        }}

        .rate-positive {{
            background: #d5f4e6;
            color: var(--production-green);
        }}

        .rate-negative {{
            background: #fadbd8;
            color: var(--loss-red);
        }}

        .rate-zero {{
            background: var(--primary-light);
            color: var(--accent-brown);
        }}

        .reactions-section {{
            margin-top: 15px;
        }}

        .reaction-item {{
            margin: 8px 0;
            padding: 15px;
            border-radius: 8px;
            font-size: 0.9em;
            border-left: 4px solid;
        }}

        .production {{
            background: linear-gradient(90deg, #e8f8f2, #f0fbf5);
            border-left-color: var(--production-green);
        }}

        .loss {{
            background: linear-gradient(90deg, #fdeaea, #fef2f2);
            border-left-color: var(--loss-red);
        }}

        .reaction-rate {{
            font-weight: 600;
            float: right;
            padding: 2px 8px;
            border-radius: 4px;
            font-size: 0.85em;
        }}

        .production .reaction-rate {{
            background: var(--production-green);
            color: white;
        }}

        .loss .reaction-rate {{
            background: var(--loss-red);
            color: white;
        }}

        .reaction-equation {{
            color: var(--accent-brown);
            font-size: 0.85em;
            margin-top: 8px;
            font-family: 'Courier New', monospace;
            background: rgba(255,255,255,0.5);
            padding: 5px 8px;
            border-radius: 4px;
        }}

        .chemical-equation {{
            color: var(--primary-dark);
            font-weight: 500;
            margin-top: 5px;
            font-size: 0.9em;
        }}

        .summary-stats {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 30px 0;
            padding: 20px;
            background: linear-gradient(135deg, var(--primary-light), #eef4ed);
            border-radius: 10px;
        }}

        .stat-item {{
            text-align: center;
        }}

        .stat-value {{
            font-size: 1.5em;
            font-weight: 600;
            color: var(--primary-dark);
        }}

        .stat-label {{
            color: var(--accent-brown);
            font-size: 0.9em;
        }}

        .toggle-button {{
            background: linear-gradient(135deg, var(--primary-medium), var(--primary-dark));
            color: white;
            border: none;
            padding: 10px 16px;
            border-radius: 6px;
            cursor: pointer;
            font-size: 0.9em;
            margin-top: 10px;
            transition: all 0.3s ease;
        }}

        .toggle-button:hover {{
            background: linear-gradient(135deg, var(--primary-dark), var(--accent-brown));
            transform: translateY(-1px);
        }}

        .reactions-details {{
            display: none;
            margin-top: 15px;
        }}

        .reactions-details.show {{
            display: block;
        }}

        .arrow {{
            color: var(--accent-orange);
            font-weight: bold;
            padding: 0 5px;
        }}

        .equilibrium-arrow {{
            color: var(--primary-medium);
            font-weight: bold;
            padding: 0 5px;
        }}
    </style>
    <script>
        function toggleReactions(speciesId) {{
            const details = document.getElementById(speciesId + '-details');
            const button = document.getElementById(speciesId + '-button');

            if (details.classList.contains('show')) {{
                details.classList.remove('show');
                button.textContent = 'Show Reactions';
            }} else {{
                details.classList.add('show');
                button.textContent = 'Hide Reactions';
            }}
        }}
    </script>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>Reaction Rates</h1>
            <div class="subtitle">Generated on {datetime.now().strftime('%B %d, %Y at %I:%M %p')}</div>
        </div>
""")

    def format_species_name(self, species):
        # Purpose: converts a species name from the internal database format to the display format.
        # Inputs are: species, a species name in the internal database format, e.g. 'HPLUS1'.
        # Outputs are: the same name converted to the display format, e.g. 'H+1'.
        # Called by: reconstruct_chemical_equation() and log_timepoint(), both in this file.
        return species.replace('PLUS', '+').replace('MINUS', '-')

    def reconstruct_chemical_equation(self, rate_law_info):
        # Purpose: turns one reaction's data back into a readable chemical equation for display.
        # Inputs are: rate_law_info, one reaction's dictionary from rate_laws (with 'consumed', 'produced', 'forward' and 'reverse' entries).
        # Outputs are: a human-readable chemical equation string, e.g. "FeII + H2O2 → FeIII + OH-", using a two-way arrow (⇌) if the reaction has both a forward and reverse rate, or an empty string if the reaction has neither reactants nor products.
        # Called by: identify_reaction_source(), in this file.
        if not rate_law_info:
            return ""

        try:
            consumed = rate_law_info.get('consumed', [])
            produced = rate_law_info.get('produced', [])

            if not consumed and not produced:
                return ""

            # Reactants are formatted for display, with whole-number coefficients shown as integers and fractional ones to 2 decimal places.
            reactants = []
            for species, coeff in consumed:
                formatted_species = self.format_species_name(species)
                if coeff == 1.0:
                    reactants.append(formatted_species)
                else:
                    if coeff == int(coeff):
                        reactants.append(f"{int(coeff)} {formatted_species}")
                    else:
                        reactants.append(f"{coeff:.2f} {formatted_species}")

            # Products are formatted the same way as reactants above.
            products = []
            for species, coeff in produced:
                formatted_species = self.format_species_name(species)
                if coeff == 1.0:
                    products.append(formatted_species)
                else:
                    if coeff == int(coeff):
                        products.append(f"{int(coeff)} {formatted_species}")
                    else:
                        products.append(f"{coeff:.2f} {formatted_species}")

            # A reaction with both a forward and a reverse rate is shown with an equilibrium arrow; one-way reactions get a plain arrow.
            forward_rate = rate_law_info.get('forward')
            reverse_rate = rate_law_info.get('reverse')

            if forward_rate and reverse_rate and reverse_rate != 'None':
                arrow = '<span class="equilibrium-arrow">⇌</span>'  # Equilibrium
            else:
                arrow = '<span class="arrow">→</span>'  # One-way

            if reactants and products:
                reactant_str = " + ".join(reactants)
                product_str = " + ".join(products)
                return f"{reactant_str} {arrow} {product_str}"
            elif reactants:  # Decomposition
                reactant_str = " + ".join(reactants)
                return f"{reactant_str} {arrow} ..."
            elif products:  # Formation from unknown
                product_str = " + ".join(products)
                return f"... {arrow} {product_str}"

        except Exception as e:
            return f"[Error reconstructing equation: {str(e)[:50]}...]"

        return ""

    def safe_pow(self, x, p):
        # Purpose: raises a number to a power without ever returning zero, NaN, or infinity, purely for display purposes.
        # Inputs are: x, a number; p, the power to raise it to.
        # Outputs are: x raised to the power p, floored at 1e-20 (rather than allowed to reach zero or become invalid), so that a rate calculated purely for display in the report never divides by zero or produces NaN/infinity.
        # Called by: log_timepoint(), in this file, when displaying each reaction's contribution to a species' rate.
        if x <= 1e-20:
            return 1e-20 if p > 0 else 1.0
        try:
            result = x ** p
            return max(result, 1e-20) if np.isfinite(result) else 1e-20
        except:
            return 1e-20

    def replace_powers_with_safe_pow(self, expr):
        # Purpose: rewrites every "**" power operation in a rate-law expression to go through safe_pow() instead.
        # Inputs are: expr, a rate-law expression as a string, e.g. containing a term like "FeII**2".
        # Outputs are: the same expression with every "**" power operation rewritten as a call to safe_pow(...), e.g. "safe_pow(FeII, 2)".
        # Called by: log_timepoint(), in this file.
        pattern = re.compile(r'([a-zA-Z_][a-zA-Z0-9_]*)\s*\*\*\s*([-+]?[0-9]*\.?[0-9]+)')  # matches a species name followed by **exponent, e.g. "FeII**2"
        def replacer(match):
            # Purpose: turns one regex match (a species and its exponent) into a safe_pow(...) call.
            # Inputs are: match, a regex match object with the species name as group 1 and the exponent as group 2.
            # Outputs are: the replacement text, e.g. "safe_pow(FeII, 2)".
            # Called by: pattern.sub(), just below, once per "**" match found in expr.
            base = match.group(1)
            exponent = match.group(2)
            return f"safe_pow({base}, {exponent})"
        return pattern.sub(replacer, expr)

    def parse_reaction_terms(self, equation_right_side):
        # Purpose: splits a species' differential equation into its individual reaction terms.
        # Inputs are: equation_right_side, the right-hand side of a species' differential equation, as a string.
        # Outputs are: a list of its top-level "+"/"-" separated terms (one per contributing reaction), ignoring any "+"/"-" that appears inside parentheses.
        # Called by: log_timepoint(), in this file.
        terms = []
        current_term = ''
        depth = 0

        for i, char in enumerate(equation_right_side):
            if char == '(':
                depth += 1
            elif char == ')':
                depth -= 1
            elif depth == 0 and char in '+-' and current_term and i > 0:
                terms.append(current_term.strip())
                current_term = char
            else:
                current_term += char

        if current_term:
            terms.append(current_term.strip())

        return terms

    def identify_reaction_source(self, term):
        # Purpose: works out which specific reaction a differential-equation term came from.
        # Inputs are: term, one term from a species' differential equation, as a string.
        # Outputs are: a 4-tuple (reaction_name, expr, rate_law_info, chemical_eq): reaction_name is a label like "Reaction 3 (forward)"; expr is the matched rate-law expression; rate_law_info is that reaction's full dictionary, or None if no match was found; chemical_eq is the human-readable equation from reconstruct_chemical_equation(), or an empty string if there was no match.
        # Called by: log_timepoint(), in this file.
        for i, rate_law in enumerate(self.rate_laws):
            forward_expr = rate_law.get('forward', '')
            reverse_expr = rate_law.get('reverse', '')

            if forward_expr and forward_expr in term:
                chemical_eq = self.reconstruct_chemical_equation(rate_law)
                return f"Reaction {i+1} (forward)", forward_expr, rate_law, chemical_eq
            elif reverse_expr and reverse_expr in term and reverse_expr != 'None':
                chemical_eq = self.reconstruct_chemical_equation(rate_law)
                return f"Reaction {i+1} (reverse)", reverse_expr, rate_law, chemical_eq

        return "Unknown reaction", term, None, ""

    def log_timepoint(self, t, y, dy_dt):
        # Purpose: writes one timepoint's section into the HTML report, including each active species' top contributing reactions.
        # Inputs are: t, the simulation time being logged; y, the concentration of every species at that time; dy_dt, their rate of change at that time.
        # Outputs are: nothing is returned; a section is appended to the HTML report for this timepoint, showing a summary and, for each active species, its top contributing reactions.
        # Called by: UniversalRobustODESystem.__call__(), in differentials.py; fenton_ode(), in main.py - both only when should_log_time() says this timepoint should be logged.
        context = {species: y[i] for i, species in enumerate(self.database_species)}
        context['safe_pow'] = self.safe_pow

        active_species = sum(1 for rate in dy_dt if abs(rate) > 1e-15)
        max_production = max(dy_dt) if len(dy_dt) > 0 else 0
        max_loss = min(dy_dt) if len(dy_dt) > 0 else 0
        total_activity = sum(abs(rate) for rate in dy_dt)

        with open(self.output_file, 'a') as f:
            f.write(f"""
        <div class="time-section">
            <div class="time-header">
                Time: {t:.2f} seconds ({t/3600:.3f} hours)
            </div>

            <div class="summary-stats">
                <div class="stat-item">
                    <div class="stat-value">{active_species}</div>
                    <div class="stat-label">Active Species</div>
                </div>
                <div class="stat-item">
                    <div class="stat-value">{max_production:.2e}</div>
                    <div class="stat-label">Max Production</div>
                </div>
                <div class="stat-item">
                    <div class="stat-value">{abs(max_loss):.2e}</div>
                    <div class="stat-label">Max Loss</div>
                </div>
                <div class="stat-item">
                    <div class="stat-value">{total_activity:.2e}</div>
                    <div class="stat-label">Total Activity</div>
                </div>
            </div>

            <div class="species-grid">
""")

            for i, species in enumerate(self.database_species):
                concentration = y[i]
                net_rate = dy_dt[i]
                formatted_name = self.format_species_name(species)
                species_id = f"species_{i}_t_{int(t)}"

                if abs(net_rate) < 1e-15:
                    rate_class = "rate-zero"
                elif net_rate > 0:
                    rate_class = "rate-positive"
                else:
                    rate_class = "rate-negative"

                f.write(f"""
                <div class="species-card">
                    <div class="species-name">{formatted_name}</div>
                    <div class="rate-info">
                        <div class="concentration">Conc: {concentration:.2e} M</div>
                        <div class="net-rate {rate_class}">{net_rate:+.2e}</div>
                    </div>
""")

                # An active species' contributing reactions are only broken down if the species is actually changing at a meaningful rate.
                if abs(net_rate) > 1e-15 and species in self.differential_equations:
                    equation = self.differential_equations[species]
                    _, right_side = equation.split('=', 1)
                    right_side = right_side.strip()

                    if right_side != '0':
                        terms = self.parse_reaction_terms(right_side)
                        term_values = []

                        for term in terms:
                            try:
                                safe_term = self.replace_powers_with_safe_pow(term)
                                value = eval(safe_term, {"safe_pow": self.safe_pow}, context)
                                term_values.append((term, value))
                            except:
                                term_values.append((term, 0.0))

                        # Terms are ranked by how much they contribute, largest first.
                        term_values.sort(key=lambda x: abs(x[1]), reverse=True)

                        # Only reactions contributing at least 1% of the species' net rate are shown, to keep the report focused on what actually matters.
                        significant_terms = [tv for tv in term_values if abs(tv[1]) > abs(net_rate) * 0.01]

                        if significant_terms:
                            f.write(f"""
                    <button class="toggle-button" id="{species_id}-button" onclick="toggleReactions('{species_id}')">
                        Show Reactions ({len(significant_terms)} contributing)
                    </button>

                    <div class="reactions-details" id="{species_id}-details">
                        <div class="reactions-section">
""")

                            for term, value in significant_terms[:10]:  # Top 10 reactions
                                reaction_type = "production" if value > 0 else "loss"
                                reaction_name, expr, rate_law_info, chemical_eq = self.identify_reaction_source(term)

                                f.write(f"""
                            <div class="reaction-item {reaction_type}">
                                <div>
                                    {reaction_name}
                                    <span class="reaction-rate">{value:+.2e}</span>
                                </div>
                                {f'<div class="chemical-equation">{chemical_eq}</div>' if chemical_eq else ''}
                                <div class="reaction-equation">{term}</div>
                            </div>
""")

                            f.write("                        </div>\n                    </div>")

                f.write("                </div>")

            f.write("""
            </div>
        </div>
""")

    def should_log_time(self, t, total_time, max_points=50):
        # Purpose: decides whether the current time is one of the handful of points that should be logged.
        # Inputs are: t, the current simulation time; total_time, the total length of the integration; max_points, how many points should be logged across the whole run (default 50).
        # Outputs are: True if t is close enough to one of max_points evenly-spaced target times (and that target hasn't already been logged), False otherwise.
        # Called by: UniversalRobustODESystem.__call__(), in differentials.py; fenton_ode(), in main.py.
        if total_time <= 0:
            return True

        target_times = np.linspace(0, total_time, max_points)

        for target_time in target_times:
            if abs(t - target_time) < total_time / (max_points * 10):  # Within 0.2% of target
                if target_time not in self.logged_times:
                    self.logged_times.add(target_time)
                    return True

        return False

    def finalise_html(self):
        # Purpose: closes off the HTML report so it's a valid, complete file.
        # Inputs are: none.
        # Outputs are: nothing is returned; the closing HTML tags are appended to the report file, and a confirmation message is printed.
        # Called by: UniversalRobustODESystem.finalise_logging(), in differentials.py; finalise_logging(), the nested function in main.py.
        with open(self.output_file, 'a') as f:
            f.write("""
    </div>
</body>
</html>""")

        print(f"Detailed reaction rates analysis saved to: {self.output_file}")
