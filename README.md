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
* fields `new_contact_rate` and `loose_contact_rate` describe 
* field `seed` aloows to fix a seed for the Pseudo-Random Number Generator (Mersenne Twister 19937) during initiation of the Contact network. Plese note that simulations performed with randomly chosen different seeds are not guaranteed to be (pseudo)independent.
* field `initial_edges` describes an initial number of edges in the Contact Network
  
Parameter `-mode` allows to run either SSATAN-X algorithm using `-SSX` or classic SSA algorithm using `-SSA`.  
  
The program is implemented for the following model: &#8594;
```S + I  I + I```
