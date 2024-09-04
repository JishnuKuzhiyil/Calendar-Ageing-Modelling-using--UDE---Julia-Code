# **Lithium-ion Battery Calendar Ageing Model using Universal Differential Equations**

This repository contains the implementation of a Lithium-ion Battery Calendar Ageing Model, developed using Universal Differential Equations (UDEs). The model is based on the article titled "*Lithium-ion Battery Degradation Modelling using Universal Differential Equations: Development of a Cost-Effective Parameterisation Methodology*" by **Jishnu Ayyangatu Kuzhiyil**, **Theodoros Damoulas**, **Ferran Brosa Planella**, and **W. Dhammika Widanage**.

## **Authors and Affiliations**

- **Jishnu Ayyangatu Kuzhiyil**  
  WMG, University of Warwick, Coventry, UK  
  Email: [jishnu-ak.ayyangatu-kuzhiyil@warwick.ac.uk](mailto:jishnu-ak.ayyangatu-kuzhiyil@warwick.ac.uk)

- **Theodoros Damoulas**  
  Department of Computer Science and Department of Statistics, University of Warwick, Coventry, UK  
  The Alan Turing Institute, London, UK

- **Ferran Brosa Planella**  
  The Faraday Institution, Quad One, Harwell Science and Innovation Campus, Didcot, UK  
  Mathematics Institute, University of Warwick, Gibbet Hill Road, Coventry, CV4 7AL, UK

- **W. Dhammika Widanage**  
  WMG, University of Warwick, Coventry, UK

## **Model Overview**

This repository provides a computational model for simulating the calendar ageing of Lithium-ion batteries. The model includes two approaches:

1. **Physics-Based Model**:  
   A Single Particle Model with Electrolyte (SPMe) forms the base electrochemical model, coupled with additional degradation differential equations to describe Solid Electrolyte Interphase (SEI) growth and associated pore clogging. The degradation is modeled using diffusion-limited and kinetically-limited SEI growth models along with a linear pore clogging model.

2. **UDE-Based Model**:  
   A Single Particle Model with Electrolyte (SPMe) forms the base electrochemical model, coupled with additional degradation differential equations are modeled as Universal Differential Equations, where neural networks are incorporated into the degradation differential equations to capture complex degradation behavior.

## Data Files

The project includes the following data files:

* `RPT_analysis_data.mat`: This file contains the capacity data used for model validation. .

* `RPTx_analysis_data.mat`: This file contains anode Loss of Active Material (LAM) data as used in the referenced article.
  
## **Usage Instructions**

To run the model, follow these steps:

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/JishnuKuzhiyil/Calendar-Ageing-Modelling-using--UDE---Julia-Code.git
   cd Calendar-Ageing-Modelling-using--UDE---Julia-Code

2. **Install the required Julia packages**:
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```

3. **Open the `Main.jl` file in your Julia environment**.

4. **Modify the following variables to simulate different calendar ageing scenarios:
   - `SOC`: Set the State of Charge (e.g., `SOC = 50` for 50% charge)
   - `Temperature`: Set the ambient temperature in degrees Celsius (e.g., `Temperature = 25`)

5. **Choose the model type by setting the `Model` variable**:
   - For the physics-based model: `Model = "Physics"`
   - For the UDE model: `Model = "UDE"`

6. **Run the `Main.jl` script**:
   ```
   julia Main.jl
   ```

## Output

Plots showing capacity and anode LAM measurements and corresponding model predictions.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

```
MIT License

Copyright (c) [2024] [Jishnu Ayyangatu Kuzhiyil]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## How to Cite

If you use this software in your research, please cite our article:

```
Kuzhiyil, J. A., et al. (2024). Lithium-ion Battery Degradation Modelling using Universal Differential Equations: Development of a Cost-Effective Parameterisation Methodology. 
```


