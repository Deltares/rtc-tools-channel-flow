# Guideline: Documenting Modelica Blocks for Automatic Read-the-Docs Generation

This project supports automatic generation of Sphinx/Read-the-Docs documentation from Modelica `Documentation(info=...)` annotations. Documentation is written once in the Modelica model and automatically converted into `.rst` pages during the documentation build.

## Overview

The documentation workflow is:

```text
Modelica (.mo)
    ↓
Documentation(info="...")
    ↓
generate_modelica_docs.py
    ↓
generated/*.rst
    ↓
Sphinx / Read the Docs
```

The Modelica file becomes the **single source of truth**.

---

# 1. Add documentation to the Modelica block

Use an `annotation(Documentation(info="..."))` block inside the model.

Example:

```modelica
model HomotopicPower
  extends Internal.PartialHomotopicPower;

  annotation (
    Documentation(info="
<html>

<p>
This model extends HomotopicVolume with a hydropower generation equation.
</p>

<p>
Power generation depends on turbine discharge and the available hydraulic head,
computed as the difference between the reservoir water level and the tailwater level.
Turbine efficiency is assumed to be constant.
</p>

<math>
P = \\eta \\rho g Q H
</math>

<p>
The optimization starts from a simplified power equation in which the hydraulic
head is assumed to be constant. This reference head must be provided as a model parameter.
/p>

</html>")
  );
end HomotopicPower;
```

---

# 2. Supported formatting

## Paragraphs

Use:

```html
<p>
This is a paragraph.
</p>
```
Generates:

```rst
This is a paragraph.
```

---

## Italics

Use:

```html
<em>Node</em>
```

Generates:

```rst
*Node*
```

---

## Bold text

Use:

```html
<strong>Important</strong>
```

Generates:

```rst
**Important**
```

---

## Lists
Use:

```html
<ul>
<li>First item</li>
<li>Second item</li>
</ul>
```
Generates:

```rst
- First item
- Second item
```

---

## Inline code

Use:

```html
<code>nout</code>
```

Generates:

```rst
``nout``
```

Useful for parameter names such as:

```html
<code>nin</code>
<code>nout</code>
<code>QIn</code>
<code>QOut</code>
```

---

## Code examples

Use multiline code blocks:

```html
<code>
Deltares.ChannelFlowSimpleRouting.Nodes.Node Alder(
  nin = 2,
  nout = 1,
  n_QForcing =0)
</code>
```

Generates:

```rst::

   Deltares.ChannelFlow.Simplerouting.Nodes.Node Alder(
     nin = 2,
     nout = 1,
     n_QForcing= 0)
```

---

## Equations

Use:
```html
<math>
h = \\frac{V}{A}
</math>
```

Generates:

```rst
.. math::

   h = \frac{V}{A}
```

The documentation generator automatically converts escaped Modelica backslashes into valid LaTeX.

---

# 3. Example: documenting a Node block

```modelica
annotation (
  Documentation(info="
<html>

<p>
A <em>node</em> connects multiple branches. It can be used to model
bifurcations and confluences in a water system.
</p>

<p>
Example:
</p>

<code>
Deltares.ChannelFlow.SimpleRouting.Nodes.Node Alder(
  nin = 2,
  nout = 1,
  nQForcing = 0,
  QIn.Q(each nominal= 0.3),
  QOut.Q(each nominal = 10)
</code>

<p>
In this example, the <em>Node</em> named <em>Alder</em> has two inflow
connectors and one outflow connector (<code>nout</code>) to represent
a confluence.
</p>
<math>
Q_{in,1} + Q_{in,2} = Q_{out}
</math>

</html>")
);
```

Generated output:

```rst
Node
====

A *Node* connects multiple branches. It can be used to model
bifurcation and confluences in a water system.

Example:

::

   Deltares.ChannelFlow.SimpleRouting.Nodes.Node Alder(
     nin = 2,
     nout = 1,
    n_QForcing = 0,
     QIn.Q(each nominal = 0.3),
     QOut.Q(each nominal = 10))

In this example, the *Node* named *Alder* has two inflow
connectors and one outflow connection (``nout``) to represent
a confluence.

.. math::

   Q_{in,1} + Q_{i,2} = Q_{out}
```

---

# 4. Package hierarchy generation

The documentation generator mirrors the Modelica package hierarchy.

Modelica:

```text
Deltares/
├── ChannelFlow/
|   ├── Hydraulic/
│   │   └── Reservoir/
│   │       ├── HomotopicPower.mo
│   │       └── HomotopicVolume.mo
│   └── SimpleRouting/│       └── Nodes/
│           └──Node.mo
```

Generated documentation:

```text
generated/
├── index.rst
├── ChannelFlow/
│   ├── index.rst
│   ├── Hydraulic/
│   │   ├── index.rst
│   │   └── Reservoir/
│   │       ├── index.rst
│   │       ├── HomotopicPower.rst
│   │       └── HomotopicVolume.rst
│   └── SimpleRouting/
│       ├── index.rst
│       └── Nodes/
│           ├── index.rst
│           └── Node.rst
```

This creates a navigation tree that matches the Modelica package structure.

---

# 5. Building the documentation

Generate the RST file:

```bash
python generate_modelica_docs.py
```

Build the documentation:

```bash
make html
```

or on windows:

```powershell
.\make.bat html
```

---

# 6. Recommended Sphinx configuration

In `conf.py`:

```python
html_theme_options = {
   "navigation_depth": 6,
}
```

This ensures the full package hierarchy remains visible in the navigation panel.

---

# 7. Recommendations

* Keep package-level documentation concise.
- Put detailed description in the Modelica model itself.
- Use `<math>` for equations.
- Use `<Node>` for parameter names and examples.
- Use `<em>` for model names and concepts.
- Avoid duplicate documentation in `.mo` and `.rst`.
- Generate documentation from Modelica whenever possible.
- Treat the Modelica file as the single source of truth.

---

# Troubleshooting

## document isn't included in any toctree"

Usually caused by:

- stale generated files;
- missing package index pages;
- missing references in generated `index.rst`.

Try deleting the generated directory and regenerating:

```bash
rm -rf source/generated
python generate_modelica_docs.py
make html
```

## Equation not rendering

Check:

```python
extensions = [
    "sphinx.ext.mathjax,
]
```

and ensure the generated RST contains:

```rst
.. math::

  h = \frac{V}{A}
```

not:

```rst
.. math::

   h = \\frac{V}{A}
```
## Inline code rendered as a code block

Use:

```html
<code>nout</code>
```

for inline code.

Use multiline `<code>` blocks only for examples:

```html
<code>
..
</code>
```
