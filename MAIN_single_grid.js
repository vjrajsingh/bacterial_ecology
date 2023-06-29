// main code file. the digestion 

// node js requirements
let Simulation = require('cacatoo') // Loads the Simulation class from local package like below
let yargs = require('yargs')
let fs = require('fs')              // for creating, editing and appending files
//const readline = require('readline'); // reads data input from the cosnole

// block for taking PG_size from console

/*const rl = readline.createInterface({
    input: process.stdin,
    output: process.stdout
});

// input public goods neighbourhood size
rl.question("Enter the neighbourhood size: ", (public_goods_neighbourhood_size) => {                                       
    if (public_goods_neighbourhood_size < 3) throw Error("Please enter a value greater than or equal to 3")
    else if (public_goods_neighbourhood_size > 11) throw Error("Error! Size of the neighbourhood should not be greater than 11")
    else if (public_goods_neighbourhood_size % 2 === 0) throw Error("Please make the size an odd number")
rl.close(); 
public_goods_neighbourhood_size = parseInt(public_goods_neighbourhood_size); // convert the string to integer*/


let rand_death = 0.1; // random death probability
let poly_loss = 0.007; // loss from polymer-digesting enzymes
let mono_loss = 0.005; // loss from monomer-digesting enzymes
let PG_size            // size of public goods neighbourhood (it  is a moore neighbourhood and should therefore be an odd number, not more than 11 and not less than 3)
let maxsteps 
let num_col = 200;
let num_row = 200;
let mix = true                     // perfect mix if true
let b_neck = false;           // implement bottleneck if true
//let PG_switch = false;        // simulates an 'operon' that checks public goods secretion under certain conditions
var aly_pioneer = 4.6;
var aly_harvester = 5.1;
var oal = 3.5;
var MM_constant_poly = 2.1;
var MM_constant_mono = 2.3;
var starting_resource_value = 1;
let data_interval = 20;
let div_threshold = 2.0;
let resource_threshold = 5.0;

// any other factors needed?

let run_num

// block for parsing arguments 

const { hideBin } = require('yargs/helpers')
const argv = yargs(hideBin(process.argv)).argv

if (typeof argv.rand_death !== "undefined") {
	rand_death = argv.rand_death
	console.log("Random Death -- user specified", rand_death)
} 
else {
	console.log("Random Death\t", rand_death)
}

if (typeof argv.poly_loss !== "undefined") {
	poly_loss = argv.poly_loss
	console.log("Loss due to aly secretion -- user specified", poly_loss)
} 
else {
	console.log("Loss due to aly secretion\t", poly_loss)
}

if (typeof argv.mono_loss !== "undefined") {
	mono_loss = argv.mono_loss
	console.log("Loss due to oal secretion -- user specified", mono_loss)
} 
else {
	console.log("Loss due to oal secretion\t", mono_loss)
}

if (typeof argv.PG_size !== "undefined") {
	PG_size = argv.PG_size
    if (argv.PG_size < 3) throw Error("Please enter a value greater than 3")
    else if (argv.PG_size > 11) throw Error("Please enter a value less than 11")
    else if (argv.PG_size % 2 === 0) throw Error("Please enter an odd number value")
	console.log("Public goods dispersal area -- user specified", PG_size)
} 
else {
    PG_size = public_goods_neighbourhood_size
	console.log("Public goods dispersal area\t", PG_size)
}

if (typeof argv.maxsteps !== "undefined") {
	maxsteps = argv.maxsteps
	console.log("Timesteps of simulation -- user specified", maxsteps)
} 
else {
	console.log("Timesteps of simulation\t", maxsteps)
}

if (typeof argv.mix !== "undefined") {
	mix = argv.mix
	console.log("Perfect Mix -- user specified", mix)
} 
else {
	console.log("Perfect Mix\t", mix)
}

if (typeof argv.b_neck !== "undefined") {
	b_neck = argv.b_neck
	console.log("Bottleneck -- user specified", b_neck)
} 
else {
	console.log("Bottleneck\t", b_neck)
}

if (typeof argv.aly_pioneer !== "undefined") {
	aly_pioneer = argv.aly_pioneer
	console.log("Aly enzyme value for pioneers -- user specified", aly_pioneer)
} 
else {
	console.log("Aly enzyme value for pioneers\t", aly_pioneer)
}

if (typeof argv.aly_harvester !== "undefined") {
	aly_harvester = argv.aly_harvester
	console.log("Aly enzyme value for harvesters -- user specified", aly_harvester)
} 
else {
	console.log("Aly enzyme value for harvesters\t", aly_harvester)
}

if (typeof argv.oal !== "undefined") {
	oal = argv.oal
	console.log("Oal enzyme value for all microbes -- user specified", oal)
} 
else {
	console.log("Oal enzyme value for all microbes\t", oal)
}

if (typeof argv.MM_constant_poly !== "undefined") {
    MM_constant_poly = argv.MM_constant_poly
    console.log("Km value for polymer digestion -- user specified", MM_constant_poly)
}
else {
    console.log("Km value for polymer digestion\t", MM_constant_poly)
}

if (typeof argv.MM_constant_mono !== "undefined") {
    MM_constant_mono = argv.MM_constant_mono
    console.log("Km value for monomer digestion -- user specified", MM_constant_mono)
}
else {
    console.log("Km value for monomer digestion\t", MM_constant_mono)
}

if (typeof argv.data_interval !== "undefined") {
	data_interval = argv.data_interval
	console.log("Timesteps between data logs -- user specified", data_interval)
} 
else {
	console.log("Timesteps between data logs\t", data_interval)
}

if (typeof argv.div_threshold !== "undefined") {
	div_threshold = argv.div_threshold
	console.log("Resources required for cells to divide -- user specified", div_threshold)
} 
else {
	console.log("Resources required for cells to divide\t", div_threshold)
}

if (typeof argv.resource_threshold !== "undefined") {
	resource_threshold = argv.resource_threshold
	console.log("Critical resource value after which resources refresh -- user specified", resource_threshold)
} 
else {
	console.log("Critical resource value after which resources refresh\t", resource_threshold)
}

if (typeof argv.run_num !== "undefined") {
    run_num = argv.run_num
    console.log("Run number -- user specified", run_num)
}
else {
    console.log("Run number\t", run_num)
}    

var dir = './data_out_'+run_num

if (!fs.existsSync(dir)){
    fs.mkdirSync(dir);
} 

let main_file = dir+"/output.dat"
if (fs.existsSync(main_file)) {
    fs.unlinkSync(main_file)
}

let run_info = dir+"/Run_Info.txt"
if(fs.existsSync(run_info)){
//	console.log("in exists block\t\t", filn)
	fs.unlinkSync(run_info)
}
fs.appendFileSync(run_info, "Run number =\t"+run_num+"\n"+"Public goods neighbourhood size =\t"+PG_size+"\n"+"Random death probability =\t"+rand_death+"\n")
fs.appendFileSync(run_info, "Maxtime =\t"+maxsteps+"\n"+"Grid size =\t"+num_col+" by "+num_row+"\n")
fs.appendFileSync(run_info, "Cost of digesting polymers =\t"+poly_loss+"\n"+"Cost of digesting monomers =\t"+mono_loss+"\n"+"Division threshold =\t"+div_threshold+"\n")
fs.appendFileSync(run_info, "Enzyme value for pioneers =\t"+aly_pioneer+"\n"+"Enzyme value for harvesters =\t"+aly_harvester+"\n"+"Enzyme value for monomer digestion (all cells) =\t"+oal+"\n")
fs.appendFileSync(run_info, "Michaelis-Menten constant (polymers) =\t"+MM_constant_poly+"\n"+"Michaelis-Menten constant (monomer) =\t"+MM_constant_mono+"\n")
fs.appendFileSync(run_info, "Perfectly mixed? =\t"+mix+"\n"+"Bottleneck? =\t"+b_neck+"\n"/*+"Public Goods Operon? =\t"+PG_switch+"\n"*/)
fs.appendFileSync(run_info, "Write out interval =\t"+data_interval+"\n"+"Resource refresh interval =\t"+resource_threshold+"\n")

//console.log(PG_size);

let config = {
    maxtime: maxsteps, // maxtime of the simulation
    ncol: num_col,
    nrow: num_row,
    seed: Math.random() * 10,
    wrap: [true, true],
    scale: 2,
    statecolours: {'type': {'p': 'red', 'h': 'yellow', 'c': 'green'}} // for display only
}

let sim = new Simulation(config)
sim.makeGridmodel('colony')       // two separate models. one for the bacteria
sim.makeGridmodel('resources')    // and one for the food particles
    
sim.initialise = function() {
    let microbes = [{type: 'p', Km_polymer: MM_constant_poly, Vmax_polymer: aly_pioneer/(PG_size * PG_size), Km_monomer: MM_constant_mono, Vmax_monomer: oal, internal_resources: starting_resource_value}, {type: 'h', Km_polymer: MM_constant_poly, Vmax_polymer: aly_harvester, Km_monomer: MM_constant_mono, Vmax_monomer: oal, internal_resources: starting_resource_value}, 
                    {type: 'c', Km_polymer: 0.0, Vmax_polymer: 0.0, Km_monomer: MM_constant_mono, Vmax_monomer: oal, internal_resources: starting_resource_value}]

 /* let pioneers = [{name:  'p', Km_polymer: 0, Vmax_polymer: 0, Km_monomer: 0, Vmax_monomer: 0, internal_resources: 50}] 
    let harvesters = [{name: 'h', Km_polymer: 0, Vmax_polymer: 0, Km_monomer: 0, Vmax_monomer: 0, internal_resources: 50}]
    let cheaters = [{name: 'c', Km_polymer: 0, Vmax_polymer: 0, Km_monomer: 0, Vmax_monomer: 0, internal_resources: 50}] */

    //sim.populateSpot(sim.colony, microbes, [0.33, 0.05, 0.34], 500, config.ncol/2, config.nrow/2)  // cannot understand the frquency array here from Bram's code
    sim.populateGrid(sim.colony, microbes, [0.33, 0.33, 0.33])
    //sim.createDisplay('colony', 'type', "")

    sim.initialGrid('resources', 'polymer_count', 100.0, 1.0) // initialise resource count. each gridpoint will have 100 units of polymer molecules
    sim.initialGrid('resources', 'monomer_count', 0.0, 1.0)   // and zero monomers
     
  //sim.createDisplay_continuous({model: 'resources', property: 'polymer_count', label: 'Polymer count', minval: 0, maxval: 100.0, fill: 'viridis'})
  //sim.createDisplay_continuous({model: 'resources', property: 'monomer_count', label: 'Monomer count', minval: 0, maxval: 1000.0, fill: 'inferno'})
}

sim.resources.refresh = function() {                                            // refreshing the old particle
    sim.initialGrid('resources', 'polymer_count', 100.0, 1.0) 
    sim.initialGrid('resources', 'monomer_count', 0.0, 1.0)
}
    
sim.colony.publicGoodsSecretion5 = function(model, col, row, direction) { // public goods neighbourhood for moore neighbourhood radius of 5
    model.moore = [[ -2, -2 ], [ -1, -2 ], [ 0, -2 ],
    [ 1, -2 ],  [ 2, -2 ],  [ -2, -1 ],
    [ -1, -1 ], [ 0, -1 ],  [ 1, -1 ],
    [ 2, -1 ],  [ -2, 0 ],  [ -1, 0 ],
    [ 0, 0 ],   [ 1, 0 ],   [ 2, 0 ],
    [ -2, 1 ],  [ -1, 1 ],  [ 0, 1 ],
    [ 1, 1 ],   [ 2, 1 ],   [ -2, 2 ],
    [ -1, 2 ],  [ 0, 2 ],   [ 1, 2 ],
    [ 2, 2 ]
    ]
    
    let x = model.moore[direction][0];
    let y = model.moore[direction][1];
    
    return model.getGridpoint(col + x, row + y);
}
    
sim.colony.publicGoodsSecretion7 = function(model, col, row, direction) { // public goods neighboirhood for moore neighbourhood radius of 7
    model.moore = [[ -3, -3 ], [ -2, -3 ], [ -1, -3 ], [ 0, -3 ],
    [ 1, -3 ],  [ 2, -3 ],  [ 3, -3 ],  [ -3, -2 ],
    [ -2, -2 ], [ -1, -2 ], [ 0, -2 ],  [ 1, -2 ],
    [ 2, -2 ],  [ 3, -2 ],  [ -3, -1 ], [ -2, -1 ],
    [ -1, -1 ], [ 0, -1 ],  [ 1, -1 ],  [ 2, -1 ],
    [ 3, -1 ],  [ -3, 0 ],  [ -2, 0 ],  [ -1, 0 ],
    [ 0, 0 ],   [ 1, 0 ],   [ 2, 0 ],   [ 3, 0 ],
    [ -3, 1 ],  [ -2, 1 ],  [ -1, 1 ],  [ 0, 1 ],
    [ 1, 1 ],   [ 2, 1 ],   [ 3, 1 ],   [ -3, 2 ],
    [ -2, 2 ],  [ -1, 2 ],  [ 0, 2 ],   [ 1, 2 ],
    [ 2, 2 ],   [ 3, 2 ],   [ -3, 3 ],  [ -2, 3 ],
    [ -1, 3 ],  [ 0, 3 ],   [ 1, 3 ],   [ 2, 3 ],
    [ 3, 3 ]
    ]
        
    let x = model.moore[direction][0];
    let y = model.moore[direction][1];
        
    return model.getGridpoint(col + x, row + y);
}
    
sim.colony.publicGoodsSecretion9 = function(model, col, row, direction) {    // public goods neighbourhood for moore neighbourhood of 9
    model.moore = [[ -4, -4 ], [ -3, -4 ], [ -2, -4 ], [ -1, -4 ], [ 0, -4 ],
    [ 1, -4 ],  [ 2, -4 ],  [ 3, -4 ],  [ 4, -4 ],  [ -4, -3 ],
    [ -3, -3 ], [ -2, -3 ], [ -1, -3 ], [ 0, -3 ],  [ 1, -3 ],
    [ 2, -3 ],  [ 3, -3 ],  [ 4, -3 ],  [ -4, -2 ], [ -3, -2 ],
    [ -2, -2 ], [ -1, -2 ], [ 0, -2 ],  [ 1, -2 ],  [ 2, -2 ],
    [ 3, -2 ],  [ 4, -2 ],  [ -4, -1 ], [ -3, -1 ], [ -2, -1 ],
    [ -1, -1 ], [ 0, -1 ],  [ 1, -1 ],  [ 2, -1 ],  [ 3, -1 ],
    [ 4, -1 ],  [ -4, 0 ],  [ -3, 0 ],  [ -2, 0 ],  [ -1, 0 ],
    [ 0, 0 ],   [ 1, 0 ],   [ 2, 0 ],   [ 3, 0 ],   [ 4, 0 ],
    [ -4, 1 ],  [ -3, 1 ],  [ -2, 1 ],  [ -1, 1 ],  [ 0, 1 ],
    [ 1, 1 ],   [ 2, 1 ],   [ 3, 1 ],   [ 4, 1 ],   [ -4, 2 ],
    [ -3, 2 ],  [ -2, 2 ],  [ -1, 2 ],  [ 0, 2 ],   [ 1, 2 ],
    [ 2, 2 ],   [ 3, 2 ],   [ 4, 2 ],   [ -4, 3 ],  [ -3, 3 ],
    [ -2, 3 ],  [ -1, 3 ],  [ 0, 3 ],   [ 1, 3 ],   [ 2, 3 ],
    [ 3, 3 ],   [ 4, 3 ],   [ -4, 4 ],  [ -3, 4 ],  [ -2, 4 ],
    [ -1, 4 ],  [ 0, 4 ],   [ 1, 4 ],   [ 2, 4 ],   [ 3, 4 ],
    [ 4, 4 ]
    ]
    
    let x = model.moore[direction][0];
    let y = model.moore[direction][1];
    
    return model.getGridpoint(col + x, row + y);
}
    
sim.colony.publicGoodsSecretion11 = function(model, col, row, direction) {  // public goods neighbourhood for moore neighbourhood of 11
    model.moore = [[-5,-5],[-4,-5],[-3,-5],[-2,-5],[-1,-5],[0,-5],[1,-5],
    [2,-5],[3,-5],[4,-5],[5,-5],[-5,-4],[-4,-4],
    [-3,-4],[-2,-4],[-1,-4],[0,-4],[1,-4],[2,-4],
    [3,-4],[4,-4],[5,-4],[-5,-3],[-4,-3],[-3,-3],
    [-2,-3],[-1,-3],[0,-3],[1,-3],[2,-3],[3,-3],
    [4,-3],[5,-3],[-5,-2],[-4,-2],[-3,-2],[-2,-2],
    [-1,-2],[0,-2],[1,-2],[2,-2],[3,-2],[4,-2],
    [5,-2],[-5,-1],[-4,-1],[-3,-1],[-2,-1],[-1,-1],
    [0,-1],[1,-1],[2,-1],[3,-1],[4,-1],[5,-1],
    [-5,0],[-4,0],[-3,0],[-2,0],[-1,0],[0,0],
    [1,0],[2,0],[3,0],[4,0],[5,0],[-5,1],
    [-4,1],[-3,1],[-2,1],[-1,1],[0,1],[1,1],
    [2,1],[3,1],[4,1],[5,1],[-5,2],[-4,2],
    [-3,2],[-2,2],[-1,2],[0,2],[1,2],[2,2],
    [3,2],[4,2],[5,2],[-5,3],[-4,3],[-3,3],
    [-2,3],[-1,3],[0,3],[1,3],[2,3],[3,3],
    [4,3],[5,3],[-5,4],[-4,4],[-3,4],[-2,4],
    [-1,4],[0,4],[1,4],[2,4],[3,4],[4,4],
    [5,4],[-5,5],[-4,5],[-3,5],[-2,5],[-1,5],
    [0,5],[1,5],[2,5],[3,5],[4,5],[5,5]
    ]
    
    let x = model.moore[direction][0];
    let y = model.moore[direction][1];
    
    return model.getGridpoint(col + x, row + y);
}

sim.colony.nextState = function(i, j) {                // next state function. this bit of code is read every time step
    let randomneigh = this.randomMoore8(this, i, j)
    let this_gp = this.grid[i][j]
    let resources_gp = sim.resources.grid[i][j]

// chance to die either randomly or due to lack of resources   
     
    if (this_gp.internal_resources < 1.0) {     // death due to lack of internal resources
        this_gp.type = 0.0
        this_gp.Km_polymer = 0.0
        this_gp.Vmax_polymer = 0.0
        this_gp.Km_monomer = 0.0
        this_gp.Vmax_monomer = 0.0
        this_gp.internal_resources = 0.0
    }

    if (this.rng.genrand_real1() < rand_death) {   // chance to die randomly
        this_gp.type = 0.0
        this_gp.Km_polymer = 0.0
        this_gp.Vmax_polymer = 0.0
        this_gp.Km_monomer = 0.0
        this_gp.Vmax_monomer = 0.0
        this_gp.internal_resources = 0.0
    }
        
    if (resources_gp.polymer_count < 0.0) {
        resources_gp.polymer_count = 0.0
    }
        
    if (resources_gp.monomer_count > 100.0) {
        resources_gp.monomer_count = 100.0
    }
    else if (resources_gp.monomer_count < 0.0) {
        resources_gp.monomer_count = 0.0
    }
    
// reproduction loop

    if (!this_gp.type && randomneigh.internal_resources > div_threshold) { // add threshold value for replication. DONE. Suitable range is 0 to 35
        this_gp.type = randomneigh.type
        this_gp.Km_polymer = randomneigh.Km_polymer
        this_gp.Vmax_polymer = randomneigh.Vmax_polymer
        this_gp.Km_monomer = randomneigh.Km_monomer
        this_gp.Vmax_monomer = randomneigh.Vmax_monomer
        randomneigh.internal_resources = this_gp.internal_resources =  randomneigh.internal_resources / 2
    }

// monomer digestion loop starts here    
        
    if (resources_gp.monomer_count > 0.0) {
        switch (this_gp.type) {
            case 'p':
                let mono_digestion_p = (this_gp.Vmax_monomer * resources_gp.monomer_count) / (this_gp.Km_monomer + resources_gp.monomer_count)
                resources_gp.monomer_count -= mono_digestion_p
                this_gp.internal_resources += mono_digestion_p
                this_gp.internal_resources -= mono_loss
            break;
            
            case 'h':
                let mono_digestion_h = (this_gp.Vmax_monomer * resources_gp.monomer_count) / (this_gp.Km_monomer + resources_gp.monomer_count)
                resources_gp.monomer_count -= mono_digestion_h
                this_gp.internal_resources += mono_digestion_h
                this_gp.internal_resources -= mono_loss
            break;

            case 'c':
                let mono_digestion_c = (this_gp.Vmax_monomer * resources_gp.monomer_count) / (this_gp.Km_monomer + resources_gp.monomer_count)
                resources_gp.monomer_count -= mono_digestion_c
                this_gp.internal_resources += mono_digestion_c
                this_gp.internal_resources -= mono_loss
            break;
        }
    }
        
// polymer digestion loop starts here

    if (resources_gp.polymer_count > 0.0 && this_gp.type === 'p') { 
        switch (PG_size) {
            case 3:                                         // public goods neighbourhood of 3
                for (let n = 0; n <= 8; n++) {
                    let poly_digestion_p = (this_gp.Vmax_polymer * this.getNeighbour(sim.resources, i, j, n).polymer_count) / (this_gp.Km_polymer + this.getNeighbour(sim.resources, i, j, n).polymer_count)
                    this.getNeighbour(sim.resources, i, j, n).polymer_count -= poly_digestion_p
                    this.getNeighbour(sim.resources, i, j, n).monomer_count += poly_digestion_p
                }
                this_gp.internal_resources -= (poly_loss * 9) // access the neighbourhood of the gridpoint and update the properties accrodingly. cost of public goods production is 9 times the initial cost   
            break;

            case 5:                                       // of 5, and so on up to 11
                for (let n = 0; n <= 24; n++) {
                    let poly_digestion_p = (this_gp.Vmax_polymer * this.publicGoodsSecretion5(sim.resources, i, j, n).polymer_count) / (this_gp.Km_polymer + this.publicGoodsSecretion5(sim.resources, i, j, n).polymer_count)
                    this.publicGoodsSecretion5(sim.resources, i, j, n).polymer_count -= poly_digestion_p
                    this.publicGoodsSecretion5(sim.resources, i, j, n).monomer_count += poly_digestion_p
                }
                this_gp.internal_resources -= (poly_loss * 25);
            break;

            case 7:
                for (let n = 0; n <= 48; n++) {
                    let poly_digestion_p = (this_gp.Vmax_polymer * this.publicGoodsSecretion7(sim.resources, i, j, n).polymer_count) / (this_gp.Km_polymer + this.publicGoodsSecretion7(sim.resources, i, j, n).polymer_count)
                    this.publicGoodsSecretion7(sim.resources, i, j, n).polymer_count -= poly_digestion_p
                    this.publicGoodsSecretion7(sim.resources, i, j, n).monomer_count += poly_digestion_p
                }
                this_gp.internal_resources -= (poly_loss * 49);
            break;
            
            case 9:
                for (let n = 0; n <= 80; n++) {
                    let poly_digestion_p = (this_gp.Vmax_polymer * this.publicGoodsSecretion9(sim.resources, i, j, n).polymer_count) / (this_gp.Km_polymer + this.publicGoodsSecretion9(sim.resources, i, j, n).polymer_count)
                    this.publicGoodsSecretion9(sim.resources, i, j, n).polymer_count -= poly_digestion_p
                    this.publicGoodsSecretion9(sim.resources, i, j, n).monomer_count += poly_digestion_p
                }
                this_gp.internal_resources -= (poly_loss * 81);
            break;
            
            case 11:
                for (let n = 0; n <= 120; n++) {
                    let poly_digestion_p = (this_gp.Vmax_polymer * this.publicGoodsSecretion11(sim.resources, i, j, n).polymer_count) / (this_gp.Km_polymer + this.publicGoodsSecretion11(sim.resources, i, j, n).polymer_count)
                    this.publicGoodsSecretion11(sim.resources, i, j, n).polymer_count -= poly_digestion_p
                    this.publicGoodsSecretion11(sim.resources, i, j, n).monomer_count += poly_digestion_p
                }
                this_gp.internal_resources -= (poly_loss * 121);
            break;    
        }
    }
        
    else if (resources_gp.polymer_count > 0.0 && this_gp.type === 'h') { // for harvester
        let poly_digestion_h = (this_gp.Vmax_polymer * resources_gp.polymer_count) / (this_gp.Km_polymer + resources_gp.polymer_count)
        resources_gp.polymer_count -= poly_digestion_h
        resources_gp.monomer_count += poly_digestion_h
        this_gp.internal_resources -= poly_loss
    }
}

sim.colony.bottleneck = function(b) {           // b is the percentage of population that gets killed by this bottleneck function
    for (let i = 0; i < this.nc; i++) {
        for (let j = 0; j < this.nr; j++) {
            if (this.rng.genrand_real1() < b) {
                // let this_gp = this.grid[i][j]
                this.grid[i][j].type = 0
                this.grid[i][j].Km_polymer = 0
                this.grid[i][j].Vmax_polymer = 0
                this.grid[i][j].Km_monomer = 0
                this.grid[i][j].Vmax_monomer = 0
                this.grid[i][j].internal_resources = 0    
            }
        }
    }
}

sim.resources.nextState = function(i, j) {} // all the update of the resources happens above in model colony, so this function is empty albeit called.

// write out functions

sim.colony.writeOut = function() {
    let time = sim.time;
    if (time % data_interval === 0) {
        let population_count = sim.colony.getPopsizes('type', ['p', 'h', 'c'])
        if (time === 0) fs.writeFileSync(main_file, ' T  '+'\t'+ ' P  '+'\t'+ ' H  '+'\t'+ ' S  ' +'\t'+'\n')
        fs.appendFileSync(main_file, time +'\t'+ population_count.join('\t') + '\n')
    }
}    

sim.colony.update = function() {
    this.asynchronous();
    this.plotPopsizes('type', ['p', 'h', 'c']); // graphs populations. only for HTML displays                         
    if (mix === true) {
        this.perfectMix();                      // perfect mixing every time step
    }

    let sum_internal_resources = 0
    for (let i = 0; i < this.nc; i++) {
        for (let j = 0; j < this.nr; j++) {
        sum_internal_resources += this.grid[i][j].internal_resources
        }
    }
    this.plotArray(["Internal resource count"], [sum_internal_resources/(this.nr*this.nc)], ["orange"], "Internal Resources") // plots the internal resource count. HTML only

    /*if (this.time % 200 === 0) {
        this.MargolusDiffusion();       // can either perform margolus diffusion consistently or every few time steps
    }*/
    
    if (b_neck === true) {
        if (this.time % 100 === 0) {
            sim.colony.bottleneck(0.99);    // calls the bottleneck function
        }
    }
}

sim.resources.update = function() {
    this.asynchronous()
    this.diffuseStates('monomer_count', 0.2)  // diffuses monomers continuously at a specified rate. NB - the rate cannot be > 0.25 as specified in Bram's code

    /*if (this.time % 100 === 0) {
        sim.initialGrid('resources', 'polymer_count', 100, 1.0)
    } */

    let sum_p = 0
    let sum_m = 0                
    for (let i = 0; i < this.nc; i++) {
        for (let j = 0; j < this.nr; j++) {
            sum_p += this.grid[i][j].polymer_count    // counts the total polymer count in the grid
            sum_m += this.grid[i][j].monomer_count    // counts the total monomer count in the grid
                
            if (sum_p < resource_threshold) {                        // introducing pristine particle whenever the last particle has critically low levels of polymer left
                sim.resources.refresh();
            }
        }                    
    }        
    this.plotArray(["Polymer count"], [sum_p/(this.nr*this.nc)], ["green"], "Polymer")
    this.plotArray(["Monomer count"], [sum_m/(this.nr*this.nc)], ['violet'], ["Monomer"])
}

sim.colony.writeOut();
sim.initialise();
sim.start();

//})
