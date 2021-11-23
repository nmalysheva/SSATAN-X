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
   
## Model
The program is implemented for the following model, described in a paper: 
* ### Contact Dynamics 
  *  Assembling of a new contact. For each pair of nodes (v<sub>i</sub> ; v<sub>i</sub>),  j &#8800; k which is not connected by an edge, the rate of assembling an edge is defined by   &lambda;<sub>j,k</sub>  &#61; &lambda;<sub>j</sub>  &lambda;<sub>k</sub> i.e. product of the assembling rates of the two nodes.
  *  Disassembling of an existing contact. For each pair of nodes (v<sub>i</sub> ; v<sub>i</sub>),  j &#8800; k that are connected by an edge, the rate of disassembling is defined as &theta;<sub>j,k</sub>  &#61; &theta;<sub>j</sub>  &theta;<sub>k</sub> i.e. the product of the disassembling rates of the two nodes.
  
* ### Epidemic dynamics
  * An infection emanating from an undiagnosed, infected individual `S` + `I` &#10230; `I` + `I` occurs with rate &gamma; > 0 if nodes j and k are connected.
  * An infection emanating from an undiagnosed, infected individual `S` + `D` &#10230; `I` + `D` occurs with rate &gamma;/2 > 0 if nodes j and k are connected.
  * An infected individual may be diagnosed with the infection `I` &#10230; `D` with rate &delta; > 0
  * Individuals may die: `I` &#10230; &#8709;, `D` &#10230; &#8709;  with rate &beta; > 0
* ### Adaptivity
  * In case of diagnosis, an individual cuts all contacts and the individual's rate of establishing new contact drops to 30% of the pre-diagnosis level, i.e &lambda;<sub>j</sub> &#61; &lambda;<sub>j</sub> &middot; 0.3. Adaprivity behaviour is inline and can be changed if desired.

