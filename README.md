# `mcp2pyscf`: Quantum Chemistry Orchestration with AI

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.x](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org/downloads/)
## 🌟 Project Overview

`mcp2pyscf` is an innovative proof-of-concept demonstrating how **Large Language Models (LLMs)**, powered by Anthropic's **Model Context Protocol (MCP)**, can intelligently orchestrate complex quantum chemistry workflows using established scientific software like **PySCF** and **RDKit**. This project aims to bridge the gap between high-level natural language requests and the intricate computational steps required for molecular simulations, paving the way for more intuitive and accessible scientific discovery.

## ✨ Features

With `mcp2pyscf`, users can interact with quantum chemistry methods through natural language prompts to perform a variety of tasks:

* **Geometry Generated for PySCF Input:** Generate ready-to-run geometry string for PySCF energy calculations with any methods from simple SMILES strings.
* **Bond Stretching Scans:** Perform and visualize energy profiles for bond stretching, generating interactive plots directly accessible in your browser.
* **Quantum Mechanical (QM) Geometry Optimization:** Conduct quantum mechanical geometry optimizations using PySCF and the `geometric` optimizer, reporting optimized energies and 3D coordinates, as well as plotting the bond stretching curve.
* **Interactive 3D Molecular Visualization:** Generate interactive 3D molecular structures in your web browser using `3Dmol.js` for immediate visual inspection.

## 🚀 Why `mcp2pyscf`? Unveiling Complexity & Impact

Computational chemistry is a powerful field, predicting molecular and material properties from first principles. However, accessing this power traditionally comes with significant hurdles. `mcp2pyscf` tackles these challenges head-on:

### **The Inherited Complexity of Electronic Structure Methods**

At its core, `mcp2pyscf` integrates with a "toolbox" of theoretical chemistry methods, each with its own mathematical underpinnings and computational demands, such as Hartree-Fock (HF), Density Functional Theory (DFT), Møller-Plesset Perturbation Theory (MP2), Coupled Cluster (CC), Configuration Interaction (CI) etc. Each of these methods requires specific inputs, careful setup, and interpretation, representing years of intricate theoretical development.

### **The Challenge of Orchestration**

Building `mcp2pyscf` involved overcoming several layers of complexity:

* **Bridging Natural Language to Code:** Translating a human request like "Optimize water's geometry" into precise, executable scientific code and tool calls.
* **Tool Interoperability (MCP):** Implementing the nascent Model Context Protocol to standardize communication between the LLM agent and diverse external Python libraries (PySCF, RDKit, Matplotlib, Geometric).
* **Scientific Library Integration:** Navigating the unique APIs, default units, and output formats of multiple specialized scientific packages. For example, ensuring `geometric` correctly interfaces with PySCF's gradient methods, and converting coordinates for `3Dmol.js` visualization.
* **Multi-Step Agentic Workflows:** Enabling the LLM to perform chained operations (e.g., "scan bond, *then* plot," or "optimize geometry, *then* visualize"), requiring sophisticated prompt engineering and robust error handling.
* **Dynamic Visualization:** Generating interactive 3D molecular structures and plots on-the-fly, delivered directly to the user's browser, bypassing the limitations of console-based AI interactions.

### **The Impact: Democratizing Quantum Chemistry**

By simplifying access to these complex computational tools, `mcp2pyscf` has the potential for significant impact across scientific and industrial domains:

* **Accelerated Materials Discovery:** Rapidly prototype and evaluate new materials for applications in batteries, catalysts, and more, significantly shortening research cycles.
* **Streamlined Drug Discovery:** Expedite the design and optimization of novel drug candidates by quickly predicting molecular properties and interactions.
* **Empowering Non-Specialists:** Lower the barrier to entry for researchers, students, and engineers to leverage advanced quantum chemical simulations without requiring deep expertise in specific software syntaxes.
* **Enhanced Research Productivity:** Automate tedious, repetitive tasks, allowing scientists to focus on higher-level analysis and innovative problem-solving.
* **The Future of Scientific Computing:** Demonstrates a paradigm shift towards intelligent, AI-orchestrated scientific workflows, setting a precedent for human-computer collaboration in research.

## 🛠️ Technical Stack

* **Large Language Model (LLM):** OpenAI's GPT-4o (or Anthropic's Claude 3.5 Sonnet, or other compatible models)
* **AI Orchestration:** Anthropic's Model Context Protocol (MCP) Python SDK
* **Quantum Chemistry:** [PySCF](https://pyscf.org/)
* **Molecular Manipulation/Geometry:** [RDKit](https://www.rdkit.org/)
* **Geometry Optimization:** [Geometric](https://geometric.readthedocs.io/en/latest/)
* **Plotting:** [Matplotlib](https://matplotlib.org/)
* **3D Visualization:** [3Dmol.js](https://3dmol.org/)
* **Core Language:** Python 3.x

## 🚀 Getting Started

To run `mcp2pyscf` locally, follow these steps:

### **Prerequisites**

* Python 3.8+
* An API key for your chosen LLM (e.g., OpenAI API key)

### **Installation**

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/](https://github.com/)[YourGitHubUsername]/mcp2pyscf.git
    cd mcp2pyscf
    ```
2.  **Create and activate a virtual environment:**
    `uv` is recommended. More info can be found in [How to Build Your Own MCP Server — Step by Step (Beginner Friendly)](https://github.com/pathintegral-institute/mcp.science/blob/main/docs/how-to-build-your-own-mcp-server-step-by-step.md)
    ```bash
    python -m venv .venv
    source .venv/bin/activate  # On macOS/Linux
    # .venv\Scripts\activate   # On Windows
    ```

### **Configuration**

1.  **Set your LLM API Key:**
    Export your API key as an environment variable:
    ```bash
    export OPENAI_API_KEY="sk-YOUR_OPENAI_KEY_HERE" # macOS/Linux
    # For Windows Command Prompt: set OPENAI_API_KEY="sk-YOUR_OPENAI_KEY_HERE"
    # For Windows PowerShell: $env:OPENAI_API_KEY="sk-YOUR_OPENAI_KEY_HERE"
    ```
    *(If using Anthropic, replace `OPENAI_API_KEY` with `ANTHROPIC_API_KEY` and update `llm_client.py` accordingly.)*

### **Usage**

Please follow [Step-by-step guide to configure MCP servers for AI client apps locally](https://github.com/pathintegral-institute/mcp.science/blob/main/docs/integrate-mcp-server-step-by-step.md) to start MCP.

2.  **Try out some prompts!**

    * `Generate a PySCF HF/STO-3G input for water.`
    * `Perform a bond stretching scan for water (O=0, H1=1) from 0.8 to 1.5 Angstroms with 10 points and plot the energy.`
    * `Perform a quantum mechanical geometry optimization for methane (C) and report the optimized energy and coordinates.`
    * `Generate a 3D geometry for ethanol (CCO) using RDKit and report the coordinates.`
    * `Visualize the water molecule (O) in 3D.`

    *(For visualization and plotting, ensure your default web browser is set up to automatically open HTML files.)*


## 🚀 Future Enhancements

* Support for more advanced quantum chemistry methods (e.g., DFT functionals, MP2, Coupled Cluster).
* Integration with external chemical databases (e.g., PubChem) for initial molecule retrieval.
* Ability to parse and interpret complex PySCF output files for richer analysis.
* More sophisticated error handling and LLM-driven debugging suggestions.
* Graphical User Interface (GUI) for a more polished user experience.

## 🤝 Contributing

Contributions are welcome! If you have ideas for improvements or encounter issues, please open a GitHub issue or submit a pull request.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgements

* A special thank you to the organizers and sponsors for [MCP x Quantum Science Hackathon](https://ai-4-science.org) for supporting this project.
* To the developers of MCP, PySCF, RDKit, Geometric, Matplotlib, and 3Dmol.js for their incredible open-source scientific software.

