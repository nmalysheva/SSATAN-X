# SSATAN-X

C++ code for implementation of the Stochastic Simulation Algorithm for effective spreading dynamics on Time-evolving Adaptive NetworX (SSATAN-X). The algorithm is described in the paper: N. Malysheva, M. von Kleist "Stochastic Simulation Algorithm for effective spreading dynamics on Time-evolving Adaptive NetworX (SSATAN-X)"  (https://biorxiv.org/cgi/content/short/2021.11.22.469498v1).    
  
The program needs LEMON (https://lemon.cs.elte.hu/trac/lemon) and JSON for Modern C++ (https://json.nlohmann.me/) Libraries to be installed.  
  
A `CMakeLists.txt` is provided for easy building.    

After compiling, the program can be called from command line using following parameters:  
```
$ SSATAN-X config.json -mode 
```
where `config.json` is a simple `json` file with settings for initial Contact network in the following format:
* field `species` describes an array of states `S`, `I`, `D` and their initial amount, as well as the death rate for the individual of given state.
* fields `new_contact_rate` and `loose_contact_rate` describe limits `[a, b]` of rates of loosing and adding a new contact. During the initialization of the 
Contact Network they are sampled from Uniform distributuin `U(a, b)`. If you want to have equal rates for each nodes, set this parameter to `[a, a]`
* field `seed` aloows to fix a seed for the Pseudo-Random Number Generator (Mersenne Twister 19937) during initiation of the Contact network. Plese note that simulations performed with randomly chosen different seeds are not guaranteed to be (pseudo)independent.
* field `initial_edges` describes an initial number of edges in the Contact Network
* field `diagnosos_rate` describes diagnosis rate in population
* field `transmission_rate` describes transmission rate in population
  
Parameter `-mode` allows to run either SSATAN-X algorithm using `-SSX` or classic SSA algorithm using `-SSA`.  
  
The program is implemented for the following model:  
Each of the nodes `i` in population assigned the rates of loosing and adding a new contact  &mdash;  &theta;<sub>i</sub>; and  &lambda;<sub>i</sub> respectively. 

|event  | rate|
| ---      | ---       |
|`S` + `I` &#10230; `I` + `I` | &gamma;|
|`S` + `D` &#10230; `I` + `D`|  &gamma;/2|
|`I` &#10230; &#8709;    |&beta;|
|`D` &#10230; &#8709; |&beta;|

Adaprivity is inline

