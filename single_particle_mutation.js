// the single-particle model with evolution of enzyme production rate and lineage tracking

// node js requirements
let Simulation = require('cacatoo') // Loads the Simulation class from local package like below
let yargs = require('yargs')
let fs = require('fs')              // for creating, editing and appending files
let MersenneTwister = require('mersenne-twister') // for seeding
let { createCanvas } = require('@napi-rs/canvas')
let { join } = require('path')
let { writeFileSync } = require('fs')

// list of variables

let death = 0.1		// particle-wide mortality rate
let birth = 0.5		// particle-wide fecundity rate
let PG_size = 3		// size of the neighbourhood for polymer digestion by broadcasters (PG_size * PG_size)
let maxsteps
let num_col = 400
let num_row = 400
let starting_biomass_value = 1.0	// biomass of the cell on initialisation
let data_interval = 100			// data inputs every few time steps
let mutation_data_interval = 80 // based on a very generalised generation time of broadcasters; to log number of mutations
let dist_interval = 10000       // time interval between entries of trait distribution
var div_threshold = 2			// biomass threshold for cell division
let resource_threshold = 10000 		// threshold for serial transfer
let starting_polymer_value = 3		// polymer value per grid point on a particle every initialisation
let diff_interval = 1			// interval between margolus diffusion
// fraction of each agent on the particle on initialisation
let broadcaster_fraction = 0.1
let tetherer_fraction = 0.1
let cheater_fraction = 0.0  // these simulations did not explicitly have cheaters
// upon exclusion of any one cell time, these variables are used
let extinction_time = 0.0
let extinction_event = false
//var enz_decay_rate = 0.5
let dT = 0.1 // [0, 1] 			// discretise time into time steps with this variable
let mt = new MersenneTwister(); 	// will create a new MersenneTwister instance for seed generation
let bottleneck_time = 0;		// to track the period between serial transfers
let bottleneck_count = 0;		
let grid_capture = false;		// boolean to trigger capturing grid images or heatmaps
let image_duration = 10;
let diffusion_rate = 0.2		// rate of oligomer diffusion
let global_mut_rate = 0.001     // probability of mutation when a daughter cell is formed
let mutation = true
let t_vmax_limit = 0.7 * dT     // limit on enzyme production rate of tetherer

let MM_polymers = 0.5 			// Km value of polymer digestion - constant for all types
let MM_monomers = 0.2 			// Km value for monomer digestion - constant for all types
let aly_b = 9 				// vmax value for broadcasters
let aly_t = 0.5 			// vmax value for tetherers
let oal_all = 0.5 			// oal value for all types
let aly_cost = 0.01 			// cost for enzyme production (proportional to enzyme produced) - constant for all types producing aly
let oal_cost = 0.01 			// cost for enzyme production (proporttional to enzyme produced) - constant for all types

// seed
// either choose your own seed or set a random one with the mersenne twister
let n = mt.random_int()

let run_num;                // replicate number

let b_to_c = 0;
let c_to_b = 0;
let t_to_c = 0;     // how many producers get mutated to cheaters due to Vmax mutating to 0
let c_to_t = 0;     // how many cheaters become broadcasters again due to Vmax mutating to greater than 0
let poly_digestion_sum = [];
let master_poly_digestion_sum = [];


// block for parsing arguments 

const { hideBin } = require('yargs/helpers')
const { start } = require('repl')
let argv = yargs(hideBin(process.argv)).argv

if (typeof argv.death !== "undefined") {
	death = argv.death
	console.log("Random Death -- user specified", death)
} 
else {
	console.log("Random Death\t", death)
}

if (typeof argv.aly_cost !== "undefined") {
	aly_cost = argv.aly_cost
	console.log("Loss due to aly secretion in broadcasters -- user specified", aly_cost)
} 
else {
	console.log("Loss due to aly secretion\t", aly_cost)
} 

if (typeof argv.oal_cost !== "undefined") {
	oal_cost = argv.oal_cost
	console.log("Loss due to oal secretion", oal_cost)
} 
else {
	console.log("Loss due to oal secretion\t", oal_cost)
}

if (typeof argv.birth !== "undefined") {
	birth = argv.birth
	console.log("Fecundity rate -- user specified", birth)
} 
else {
	console.log("Fecundity rate\t", birth)
}

if (typeof argv.PG_size !== "undefined") {
	PG_size = argv.PG_size
    if (argv.PG_size < 3) throw Error("Please enter a value greater than 3")
    else if (argv.PG_size > 11) throw Error("Please enter a value less than 11")
    else if (argv.PG_size % 2 === 0) throw Error("Please enter an odd number value")
	console.log("Public goods dispersal area -- user specified", PG_size)
} 
else {
	console.log("Public goods dispersal area\t", PG_size)
}

if (typeof argv.maxsteps !== "undefined") {
	maxsteps = argv.maxsteps
	console.log("Timesteps of simulation -- user specified", maxsteps)
} 
else {
	console.log("Timesteps of simulation\t", maxsteps)
}

if (typeof argv.broadcaster_fraction !== "undefined") {
	broadcaster_fraction = argv.broadcaster_fraction
	console.log("Initial broadcaster fraction -- user specified", broadcaster_fraction)
} 
else {
	console.log("Initial broadcaster fraction\t", broadcaster_fraction)
}

if (typeof argv.tetherer_fraction !== "undefined") {
	tetherer_fraction = argv.cheater_fraction
	console.log("Initial tetherer fraction -- user specified", tetherer_fraction)
} 
else {
	console.log("Initial tetherer fraction\t", tetherer_fraction)
}

if (typeof argv.aly_b !== "undefined") {
	aly_b = argv.aly_b
	console.log("Aly enzyme value for broadcasters -- user specified", aly_b)
} 
else {
	console.log("Aly enzyme value for broadcasters\t", aly_b)
}

if (typeof argv.aly_t !== "undefined") {
	aly_t = argv.aly_t
	console.log("Aly enzyme value for tetherers -- user specified", aly_t)
} 
else {
	console.log("Aly enzyme value for tetherers\t", aly_t)
}

if (typeof argv.oal_all !== "undefined") {
	oal_all = argv.oal_all
	console.log("Oal enzyme value for all microbes -- user specified", oal_all)
} 
else {
	console.log("Oal enzyme value for all microbes\t", oal_all)
}

if (typeof argv.MM_polymers !== "undefined") {
    MM_polymers = argv.MM_polymers
    console.log("Km value for polymer digestion in broadcasters -- user specified", MM_polymers)
}
else {
    console.log("Km value for polymer digestion\t", MM_polymers)
}

if (typeof argv.MM_monomers !== "undefined") {
    MM_monomers = argv.MM_monomers
    console.log("Km value for monomers digestion -- user specified", MM_monomers)
}
else {
    console.log("Km value for monomer digestion\t", MM_monomers)
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

if (typeof argv.starting_polymer_value !== "undefined") {
	starting_polymer_value = argv.starting_polymer_value
	console.log("Initial polymer count per grid point -- user specified", starting_polymer_value)
} 
else {
	console.log("Initial polymer count per grid point\t", starting_polymer_value)
}

if (typeof argv.diffusion_rate !== "undefined") {
	diffusion_rate = argv.diffusion_rate
	console.log("Oligomer diffusion rate -- user specified", diffusion_rate)
} 
else {
	console.log("Oligomer diffusion rate\t", diffusion_rate)
}

if (typeof argv.diff_interval !== "undefined") {
	diff_interval = argv.diff_interval
	console.log("Cell mobility -- user specified", diff_interval)
} 
else {
	console.log("Cell mobility\t", diff_interval)
}

if (typeof argv.global_mut_rate !== "undefined") {
	global_mut_rate = argv.global_mut_rate
	console.log("Global mutation rate -- user specified", global_mut_rate)
} 
else {
	console.log("Global mutation rate\t", global_mut_rate)
}

if (typeof argv.run_num !== "undefined") {
    run_num = argv.run_num
    console.log("Run number -- user specified", run_num)
}
else {
    console.log("Run number\t", run_num)
}

if (typeof argv.dT !== "undefined") {
    dT = argv.dT
    console.log("dT value -- user specified", dT)
}
else {
    console.log("dT value\t", dT)
}

// all of my boolean variables I want passed as comamnd line arguments (without any JS bullshit of type coercion) are handled here

argv = require("yargs")
.options({
    'mar_diffusion': {
      describe: 'Specify mar_diffusion',
      type: 'boolean',
      default: true,
    },
    'mix': {
      describe: 'Specify mix',
      type: 'boolean',
      default: false,
    },
    'b_neck': {
      describe: 'Specify b_neck',
      type: 'boolean',
      default: true,
    },
    'take_txt': {
      describe: 'Flag to specify whether to take txt file',
      type: 'boolean',
      default: true,
    },
    'take_png': {
      describe: 'Flag to specify whether to take png file',
      type: 'boolean',
      default: false,
    },
  })
  .argv;

// now these options can be accessed as properties of the argv object
let mar_diffusion = argv.mar_diffusion;
let mix = argv.mix;
let b_neck = argv.b_neck;
let take_txt = argv.take_txt;
let take_png = argv.take_png;

console.log('Margolus Diffusion?:', mar_diffusion);
console.log('Perfectly Mixed?:', mix);
console.log('Bottleneck?:', b_neck);
console.log('Making text files of the grid?:', take_txt);
console.log('Making png files of the grid in runtime?:', take_png);

let home = 'output/p_'+aly_b+'_h_'+aly_t+'_mobility_'+diff_interval+'_diff_'+diffusion_rate         // example address name. change this according to your needs. i preferred to name directories based on parameters being tested
if (!fs.existsSync(home)){
    fs.mkdirSync(home);
}

var dir = home+'/data_out_'+run_num

if (!fs.existsSync(dir)){
    fs.mkdirSync(dir);
} 

let main_file = dir+"/output.dat"
if (fs.existsSync(main_file)) {
    fs.unlinkSync(main_file)
}

let Vmax_file = dir+"/Vmax_file.dat"
if (fs.existsSync(Vmax_file)) {
    fs.unlinkSync(Vmax_file)
}

let heatmaps = dir+"/heatmaps"
if (!fs.existsSync(heatmaps)) {
    fs.mkdirSync(heatmaps)
}

let movie = dir+"/images"
if (!fs.existsSync(movie)) {
    fs.mkdirSync(movie)
}

let Vmax_dist_b = dir+"/vmax_distribution_b.csv"
if (fs.existsSync(Vmax_dist_b)) {
    fs.unlinkSync(Vmax_dist_b)
}

let Vmax_dist_t = dir+"/vmax_distribution_t.csv"
if (fs.existsSync(Vmax_dist_t)) {
    fs.unlinkSync(Vmax_dist_t)
}

let lineage_dist_b = dir+"/lineage_distribution_b.csv"
if (fs.existsSync(lineage_dist_b)) {
    fs.unlinkSync(lineage_dist_b)
}

let lineage_dist_t = dir+"/lineage_distribution_t.csv"
if (fs.existsSync(lineage_dist_t)) {
    fs.unlinkSync(lineage_dist_t)
}

let bias_file = dir+"/bias_check.dat"
if (fs.existsSync(bias_file)) {
    fs.unlinkSync(bias_file)
}

let resource_count = dir+"/resource_count.dat"
if (fs.existsSync(resource_count)) {
    fs.unlinkSync(resource_count)
}

let bottleneck_time_file = dir+"/bottleneck_time.dat"
if (fs.existsSync(bottleneck_time_file)) {
    fs.unlinkSync(bottleneck_time_file)
}

let monomer_count = dir+"/monomer_count.dat"
if (fs.existsSync(monomer_count)) {
    fs.unlinkSync(monomer_count)
}

let run_info = dir+"/Run_Info.txt"
if(fs.existsSync(run_info)){
//	console.log("in exists block\t\t", filn)
	fs.unlinkSync(run_info)
}

let aly_broadcaster = aly_b * dT;
let aly_tetherer = aly_t * dT;
let oal = oal_all * dT;
let rand_death = death * dT;
let f = birth * dT;

// write all the simulation information (parameters used, variables tested, etc into a text file)

fs.appendFileSync(run_info, "All values are ignoring dT"+"\n")
fs.appendFileSync(run_info, "Run number =\t"+run_num+"\n"+"Public goods neighbourhood size =\t"+PG_size+"\n"+"Random death probability =\t"+death+"\n"+"Fecundity rate =\t"+birth+"\n")
fs.appendFileSync(run_info, "Maxtime =\t"+maxsteps+"\n"+"Grid size =\t"+num_col+" by "+num_row+"\n")
fs.appendFileSync(run_info, "Cost of digesting polymers in broadcasters =\t"+aly_cost+"\n"+"Cost of digesting monomers in broadcasters =\t"+oal_cost+"\n"+"Division threshold =\t"+div_threshold+"\n")
fs.appendFileSync(run_info, "Enzyme value for broadcasters =\t"+aly_b+"\n"+"Enzyme value for tetherers =\t"+aly_t+"\n"+"Enzyme value for monomer digestion (all cells) =\t"+oal_all+"\n")
fs.appendFileSync(run_info, "Michaelis-Menten constant (polymer) =\t"+MM_polymers+"\n"+"Michaelis-Menten constant (monomer) =\t"+MM_monomers+"\n")
fs.appendFileSync(run_info, "Perfectly mixed? =\t"+mix+"\n"+"Bottleneck? =\t"+b_neck+"\n"+"Margolus Diffusion? =\t"+mar_diffusion+"\n"+"Cell mobility =\t"+diff_interval+"\n"/*+"Public Goods Operon? =\t"+PG_switch+"\n"*/)
fs.appendFileSync(run_info, "Write out interval =\t"+data_interval+"\n"+"Resource refresh interval =\t"+resource_threshold+"\n"+"Polymer units per grid point =\t"+starting_polymer_value+"\n"+"Oligomer diffusion rate =\t"+diffusion_rate+"\n")
fs.appendFileSync(run_info, "Ratio of broadcaster to tetherer =\t"+broadcaster_fraction+":"+tetherer_fraction+"\n")
fs.appendFileSync(run_info, "Taking PNG?"+"\t"+take_png+"\n"+"Taking TXT?"+"\t"+take_txt+"\n")

let config = {
    maxtime: maxsteps,
    ncol: num_col,
    nrow: num_row,
    seed: n,
    wrap: [true, true],
    scale: 3,
    statecolours: {'type':{'b': 'amaranth', 't': 'gold', 'c': 'cyan'}}
}

fs.appendFileSync(run_info, "seed =\t"+config.seed+"\n")
let sim = new Simulation(config)
sim.makeGridmodel('colonies')
sim.makeGridmodel('resources')
  
let microbes = [{type: 'b', Vmax_polymer: aly_broadcaster/Math.pow(PG_size, 2), Km_polymer: MM_polymers, Vmax_monomer: oal, Km_monomer: MM_monomers, cell_biomass: starting_biomass_value, uID: 0},
	            {type: 't', Vmax_polymer: aly_tetherer, Km_polymer: MM_polymers, Vmax_monomer: oal, Km_monomer: MM_monomers, cell_biomass: starting_biomass_value, uID: 0},
                {type: 'c', Vmax_polymer: 0, Km_polymer: 0, Vmax_monomer: oal, Km_monomer: MM_monomers, cell_biomass: starting_biomass_value, uID: 0}]
  
sim.initialise = function() {
    sim.populateGrid(sim.colonies, microbes, [broadcaster_fraction, tetherer_fraction, cheater_fraction])       // initialise the particle with cell types
    //sim.createDisplay('colonies', 'type', "")
      
    sim.initialGrid('resources', 'polymer_count', starting_polymer_value, 1.0)                                  // initialise the particle with resources
    sim.initialGrid('resources', 'monomer_count', 0.0, 1.0)
      
    //sim.createDisplay_continuous({model: 'resources', property: 'polymer_count', label: 'Polymer count', minval: 0, maxval: 3.0, fill: 'viridis'})
    //sim.createDisplay_continuous({model: 'resources', property: 'monomer_count', label: 'Monomer count', minval: 0, maxval: 3.0, fill: 'inferno'})
}
 
sim.resources.refresh = function() {                                                            // function for creating a fresh, uncolonised particle with no oligomers
    sim.initialGrid('resources', 'polymer_count', starting_polymer_value, 1.0)
    sim.initialGrid('resources', 'monomer_count', 0.0, 1.0)
}

// using node's canvas module to capture grid images (optional)
const canvasWidth = num_row * config.scale;
const canvasHeight = num_col * config.scale;
const cellSize = config.scale;
let canvas = createCanvas(canvasWidth, canvasHeight);
let ctx = canvas.getContext('2d');
let buf;
let color;

// additional public goods neighbourhood sizes available

sim.colonies.publicGoodsSecretion5 = function(model, col, row, direction) { // public goods neighbourhood for moore neighbourhood radius of 5
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
    
sim.colonies.publicGoodsSecretion7 = function(model, col, row, direction) { // public goods neighboirhood for moore neighbourhood radius of 7
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
    
sim.colonies.publicGoodsSecretion9 = function(model, col, row, direction) {    // public goods neighbourhood for moore neighbourhood of 9
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
    
sim.colonies.publicGoodsSecretion11 = function(model, col, row, direction) {  // public goods neighbourhood for moore neighbourhood of 11
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

// the nextState function iterates over every grid point and updates them accordingly
 
sim.colonies.nextState = function(i, j) {
    let randomneigh = this.randomMoore8(this, i, j)
    let this_gp = this.grid[i][j]
    let resources_gp = sim.resources.grid[i][j]
    let poly_digestion_total = 0;

// chance to die either randomly or due to lack of resources   
     
    if (this_gp.cell_biomass < 0.0) {     // death due to lack of biomass
        this_gp.type = 0.0
        this_gp.Km_polymer = 0.0
        this_gp.Vmax_polymer = 0.0
        this_gp.Km_monomer = 0.0
        this_gp.Vmax_monomer = 0.0
        this_gp.cell_biomass = 0.0
    }
    
    if (sim.rng.genrand_real1() < rand_death) {     // death due to mortality
        this_gp.type = 0.0
        this_gp.Km_polymer = 0.0
        this_gp.Vmax_polymer = 0.0
        this_gp.Km_monomer = 0.0
        this_gp.Vmax_monomer = 0.0
        this_gp.cell_biomass = 0.0
    }
    
// reproduction loop

    if (mutation === true) {
        if (!this_gp.type && randomneigh.internal_resources > div_threshold) {
            if (this.rng.genrand_real1() < f) {
                this_gp.type = randomneigh.type
                this_gp.Km_polymer = randomneigh.Km_polymer
                this_gp.Km_monomer = randomneigh.Km_monomer
                this_gp.Vmax_monomer = randomneigh.Vmax_monomer
                this_gp.Vmax_polymer = randomneigh.Vmax_polymer
                this_gp.uID = randomneigh.uID
                randomneigh.internal_resources = this_gp.internal_resources = randomneigh.internal_resources / 2

                if (this.rng.genrand_real1() < global_mut_rate) {       // check if mutation will occur
                    //console.log("mutation")   
                    let m_amount = (this.rng.genrand_real1() - 0.5) * 1.0 * dT; // generates a random number between -0.5 and 0.5
                    //console.log(m_amount)
                    this_gp.Vmax_polymer += m_amount
                    this_gp.Vmax_polymer = Math.max(this_gp.Vmax_polymer, 0)
                    if (this_gp.type === 't' && this_gp.Vmax_polymer > t_vmax_limit) {this_gp.Vmax_polymer = t_vmax_limit}
                    if (randomneigh.Vmax_polymer === 0 && randomneigh.type === 'b' && this_gp.Vmax_polymer > randomneigh.Vmax_polymer) {    // restrict vmax of tetherers to the set limit
                        c_to_b++
                        //console.log(c_to_b)
                    }
                    else if (randomneigh.Vmax_polymer === 0 && randomneigh.type === 't' && this_gp.Vmax_polymer > randomneigh.Vmax_polymer) {
                        c_to_t++;
                    }
                    else if (randomneigh.Vmax_polymer > 0 && randomneigh.type === 'b' && this_gp.Vmax_polymer === 0) {
                        b_to_c++;
                    }        
                    else if (randomneigh.Vmax_polymer > 0 && randomneigh.type === 't' && this_gp.Vmax_polymer === 0) {
                        t_to_c++
                        //console.log(b_to_c)
                    }
                }
            }
        }
    }

    else {
        if (!this_gp.type && randomneigh.cell_biomass > div_threshold) {
            if (this.rng.genrand_real1() < f) {
                this_gp.type = randomneigh.type
                this_gp.Km_polymer = randomneigh.Km_polymer
                this_gp.Vmax_polymer = randomneigh.Vmax_polymer
                this_gp.Km_monomer = randomneigh.Km_monomer
                this_gp.Vmax_monomer = randomneigh.Vmax_monomer
                randomneigh.cell_biomass = this_gp.cell_biomass =  randomneigh.cell_biomass / 2
            }
        }    
    }

// monomer digestion loop starts here    
        
    if (resources_gp.monomer_count > 0.0) {
		switch (this_gp.type) {
            case 'c':
                let mono_digestion_c = (this_gp.Vmax_monomer * resources_gp.monomer_count) / (this_gp.Km_monomer + resources_gp.monomer_count)
                if (resources_gp.monomer_count < mono_digestion_c) {resources_gp.monomer_count = mono_digestion_c}
                resources_gp.monomer_count -= mono_digestion_c
                this_gp.cell_biomass += mono_digestion_c
                this_gp.cell_biomass -= oal_cost * this_gp.Vmax_monomer
            break;

            case 't':
                let mono_digestion_t = (this_gp.Vmax_monomer * resources_gp.monomer_count) / (this_gp.Km_monomer + resources_gp.monomer_count)
                if (resources_gp.monomer_count < mono_digestion_t) {resources_gp.monomer_count = mono_digestion_t}
                resources_gp.monomer_count -= mono_digestion_t
                this_gp.cell_biomass += mono_digestion_t
                this_gp.cell_biomass -= oal_cost * this_gp.Vmax_monomer
            break;

			case 'b':
				let mono_digestion_b = (this_gp.Vmax_monomer * resources_gp.monomer_count) / (this_gp.Km_monomer + resources_gp.monomer_count)
                if (resources_gp.monomer_count < mono_digestion_b) {resources_gp.monomer_count = mono_digestion_b}
                resources_gp.monomer_count -= mono_digestion_b
                this_gp.cell_biomass += mono_digestion_b
                this_gp.cell_biomass -= oal_cost * this_gp.Vmax_monomer
			break;	
        }
    }
        
// polymer digestion loop starts here

    if (resources_gp.polymer_count > 0.0 && this_gp.type === 't') {
        let poly_digestion_t = (this_gp.Vmax_polymer * resources_gp.polymer_count) / (this_gp.Km_polymer + resources_gp.polymer_count)
        if (resources_gp.polymer_count < poly_digestion_t) {resources_gp.polymer_count = poly_digestion_t}
        resources_gp.polymer_count -= poly_digestion_t
        resources_gp.monomer_count += poly_digestion_t
        poly_digestion_total += poly_digestion_t
        this_gp.cell_biomass -= aly_cost * this_gp.Vmax_polymer
        poly_digestion_sum.push(poly_digestion_total)
    }

    else if (resources_gp.polymer_count > 0.0 && this_gp.type === 'b') {
        let poly_digestion_b = []
        switch (PG_size) {
            case 3:                                         // public goods neighbourhood of 3
                for (let n = 0; n <= 8; n++) {
                    poly_digestion_b[n] = (this_gp.Vmax_polymer * this.getNeighbour(sim.resources, i, j, n).polymer_count) / (this_gp.Km_polymer + this.getNeighbour(sim.resources, i, j, n).polymer_count)
                    if (this.getNeighbour(sim.resources, i, j, n).polymer_count < poly_digestion_b[n]) {this.getNeighbour(sim.resources, i, j, n).polymer_count = poly_digestion_b[n]}
                    this.getNeighbour(sim.resources, i, j, n).polymer_count -= poly_digestion_b[n]    // access the neighbourhood of the gridpoint and update the properties accrodingly. 
                    this.getNeighbour(sim.resources, i, j, n).monomer_count += poly_digestion_b[n]
                    poly_digestion_total += poly_digestion_b[n]
                }
                //console.log(poly_digestion_total)
                poly_digestion_sum.push(poly_digestion_total)
                //console.log(poly_digestion_p)
                //console.log(poly_digestion_p)
                this_gp.cell_biomass -= (aly_cost * this_gp.Vmax_polymer * Math.pow(PG_size, 2)) // cost of public goods production is constant cost times the enzyme units secreted.
            break; 
        
            case 5:                                       // of 5, and so on up to 11
                for (let n = 0; n <= 24; n++) {
                    poly_digestion_b[n] = (this_gp.Vmax_polymer * this.publicGoodsSecretion5(sim.resources, i, j, n).polymer_count) / (this_gp.Km_polymer + this.publicGoodsSecretion5(sim.resources, i, j, n).polymer_count)
                    if (this.publicGoodsSecretion5(sim.resources, i, j, n).polymer_count < poly_digestion_b[n]) {this.publicGoodsSecretion5(sim.resources, i, j, n).polymer_count = poly_digestion_b[n]}
                    this.publicGoodsSecretion5(sim.resources, i, j, n).polymer_count -= poly_digestion_b[n]
                    this.publicGoodsSecretion5(sim.resources, i, j, n).monomer_count += poly_digestion_b[n]
                }
                    
                this_gp.cell_biomass -= (aly_cost * this_gp.Vmax_polymer * Math.pow(PG_size, 2));
            break;
        
            case 7:
                for (let n = 0; n <= 48; n++) {
                    poly_digestion_b[n] = (this_gp.Vmax_polymer * this.publicGoodsSecretion7(sim.resources, i, j, n).polymer_count) / (this_gp.Km_polymer + this.publicGoodsSecretion7(sim.resources, i, j, n).polymer_count)
                    if (this.publicGoodsSecretion7(sim.resources, i, j, n).polymer_count < poly_digestion_b[n]) {this.publicGoodsSecretion7(sim.resources, i, j, n).polymer_count = poly_digestion_b[n]}
                    this.publicGoodsSecretion7(sim.resources, i, j, n).polymer_count -= poly_digestion_b[n]
                    this.publicGoodsSecretion7(sim.resources, i, j, n).monomer_count += poly_digestion_b[n]   
                }
                    
                this_gp.cell_biomass -= (aly_cost * this_gp.Vmax_polymer * Math.pow(PG_size, 2));
            break;
                    
            case 9:
                for (let n = 0; n <= 80; n++) {
                    poly_digestion_b[n] = (this_gp.Vmax_polymer * this.publicGoodsSecretion9(sim.resources, i, j, n).polymer_count) / (this_gp.Km_polymer + this.publicGoodsSecretion9(sim.resources, i, j, n).polymer_count)
                    if (this.publicGoodsSecretion9(sim.resources, i, j, n).polymer_count < poly_digestion_b[n]) {this.publicGoodsSecretion9(sim.resources, i, j, n).polymer_count = poly_digestion_b[n]}
                    this.publicGoodsSecretion9(sim.resources, i, j, n).polymer_count -= poly_digestion_b[n]
                    this.publicGoodsSecretion9(sim.resources, i, j, n).monomer_count += poly_digestion_b[n]
                }
                    
                this_gp.cell_biomass -= (aly_cost * this_gp.Vmax_polymer * Math.pow(PG_size, 2));
            break;
                    
            case 11:
                for (let n = 0; n <= 120; n++) {
                    poly_digestion_b[n] = (this_gp.Vmax_polymer * this.publicGoodsSecretion11(sim.resources, i, j, n).polymer_count) / (this_gp.Km_polymer + this.publicGoodsSecretion11(sim.resources, i, j, n).polymer_count)
                    if (this.publicGoodsSecretion11(sim.resources, i, j, n).polymer_count < poly_digestion_b[n]) {this.publicGoodsSecretion11(sim.resources, i, j, n).polymer_count = poly_digestion_b[n]}
                    this.publicGoodsSecretion11(sim.resources, i, j, n).polymer_count -= poly_digestion_b[n]
                    this.publicGoodsSecretion11(sim.resources, i, j, n).monomer_count += poly_digestion_b[n]
                }
                    
                this_gp.cell_biomass -= (aly_cost * this_gp.Vmax_polymer * Math.pow(PG_size, 2));
            break;
        }
    }
}


sim.colonies.bottleneck = function(f) {           // f is the fraction of the community that gets removed by this bottleneck function
    for (let i = 0; i < this.nc; i++) {
        for (let j = 0; j < this.nr; j++) {
            if (this.grid[i][j].type !== 0) {
                if (this.rng.genrand_real1() < f) {
                    // let this_gp = this.grid[i][j]
                    this.grid[i][j].type = 0
                    this.grid[i][j].Km_polymer = 0
                    this.grid[i][j].Vmax_polymer = 0
                    this.grid[i][j].Km_monomer = 0
                    this.grid[i][j].Vmax_monomer = 0
                    this.grid[i][j].cell_biomass = 0
                }
            }
        }
    }
}

sim.resources.nextState = function(i, j) {}

sim.colonies.writeOut = function() {
    let time = sim.time;
    if (time % data_interval === 0) {
    let population_count = sim.colonies.getPopsizes('type', ['b', 't', 'c'])
    if (this.time === 0) fs.writeFileSync(main_file, ' T  '+'\t'+ ' B  '+'\t'+ ' T  '+'\t'+ ' C  '+'\n')
    fs.appendFileSync(main_file, time +'\t'+ population_count.join('\t') + '\n')
    }
}

// function to generate heatmaps for the grid with biomass and cell type

sim.colonies.heatmap_type_biomass = function() {
    let heatmap_file = heatmaps+"/heatmap_"+this.time+".txt"	
    for (let i = 0; i < sim.colonies.nc; i++) {
        for (let j = 0; j < sim.colonies.nr; j++) {
            if (sim.colonies.grid[i][j].type !== 0) {
                let point = i+'\t'+j+'\t'+sim.colonies.grid[i][j].type+'\t'+sim.colonies.grid[i][j].cell_biomass.toExponential(2)+'\n'
                fs.appendFileSync(heatmap_file,point)
            }
            else if (!sim.colonies.grid[i][j].type) {
                let point = i+'\t'+j+'\t'+sim.colonies.grid[i][j].type+'\t'+sim.colonies.grid[i][j].cell_biomass.toExponential(2)+'\n'
                fs.appendFileSync(heatmap_file,point)
            }	
        }
    }	
}

/*

one can also use the heatmap code block to generate txt files reflecting the numerical enzyme production rate values on the grid points
however, remember to do the appropriate operations concerning the dT variable

*/

sim.resources.heatmap_polymer = function() {
    let heatmap_file = heatmaps+"/heatmap_polymer_"+this.time+".txt"	
    for (let i = 0; i < this.nc; i++) {
        for (let j = 0; j < this.nr; j++) {
            let point = i+'\t'+j+'\t'+this.grid[i][j].polymer_count.toExponential(2)+'\n'
            fs.appendFileSync(heatmap_file,point)	
        }
    }	
}

sim.resources.heatmap_monomer = function() {
    let heatmap_file = heatmaps+"/heatmap_monomer_"+this.time+".txt"
    for (let i = 0; i < this.nc; i++) {
        for (let j = 0; j < this.nr; j++) {
            let point = i+'\t'+j+'\t'+this.grid[i][j].monomer_count.toExponential(2)+'\n'
            fs.appendFileSync(heatmap_file,point)
        }
    }
}

sim.colonies.update = function() {
    this.asynchronous();
    if (sim.time % (diff_interval) === 0 && (mar_diffusion === true)) {
        (this.MargolusDiffusion())
    }
    else if (mix === true) {
        (this.perfectMix())
    }

    if (this.time === 100000) {                         // at this time step, assign uIDs based on current trait values
        for (let i = 0; i < this.nc; i++) {
            for (let j = 0; j < this.nr; j++) {
                if (this.grid[i][j].type === 'b') {
                    if (this.grid[i][j].Vmax_polymer === 0.0) {
                        this.grid[i][j].uID = 'b_c';
                        
                    } else {
                        let lowerRange = Math.floor(this.grid[i][j].Vmax_polymer / dT * Math.pow(PG_size, 2));
                        let upperRange = lowerRange + 1;
                        this.grid[i][j].uID = `${lowerRange} to ${upperRange}`;
                        //console.log(this.grid[i][j].uID)
                    }
                }
                else if (this.grid[i][j].type === 't') {
                        if (this.grid[i][j].Vmax_polymer === 0.0) {
                            this.grid[i][j].uID = 't_c';
                            
                        } else {
                            let lowerRange = Math.floor(this.grid[i][j].Vmax_polymer / dT);
                            let upperRange = lowerRange + 1;
                            this.grid[i][j].uID = `${lowerRange} to ${upperRange}`;
                            //console.log(this.grid[i][j].uID)
                        }
                }
            }        
        }
    }

    let sum_b = 0
    let sum_cell_biomass_b = 0
    let avg_biomass_b = 0
	let sum_t = 0
	let sum_cell_biomass_t = 0
	let avg_biomass_t = 0
    let sum_c = 0
    let sum_cell_biomass_c = 0
    let avg_biomass_c = 0
    let sum_Vmax_b = 0
    let sum_Vmax_t = 0
    let avg_Vmax_b = 0
    let avg_Vmax_t = 0

    let vmax_dist_b = []
    let vmax_value_b = 0.0
    let vmax_dist_t = []
    let vmax_value_t = 0.0
    let uid_distribution_b = {};
    let uid_distribution_t = {};
    
    for (let i = 0; i < this.nc; i++) {
        for (let j = 0; j < this.nr; j++) {
            if (this.grid[i][j].type === 'b') {
                sum_b++;
                sum_cell_biomass_b += this.grid[i][j].cell_biomass
                sum_Vmax_b += this.grid[i][j].Vmax_polymer
            }
			else if (this.grid[i][j].type === 't') {
                sum_t++;
                sum_cell_biomass_t += this.grid[i][j].cell_biomass
                sum_Vmax_t += this.grid[i][j].Vmax_polymer
            }
            else if (this.grid[i][j].type === 'c') {
                sum_c++;
                sum_cell_biomass_c += this.grid[i][j].cell_biomass
            }
        }          
    }

    avg_biomass_b = sum_cell_biomass_b/sum_b
    avg_biomass_t = sum_cell_biomass_t/sum_t
    avg_biomass_c = sum_cell_biomass_c/sum_c
    avg_Vmax_b = (sum_Vmax_b/sum_b)/dT * Math.pow(PG_size, 2)
    avg_Vmax_t = (sum_Vmax_t/sum_t)/dT

    //this.plotArray(["Internal resources in broadcasters"], [avg_biomass_p], ["red"], "Internal resources in broadcasters") // plots the internal resource count in broadcasters. HTML only
	//this.plotArray(["Internal resources in tetherers"], [avg_biomass_h], ["yellow"], "Internal resources in tetherers") // plots the internal resource count in tetherers. HTML only

    sim.colonies.writeOut()
    if (this.time % data_interval === 0) {
		let out_1 = (avg_biomass_b.toExponential(2)) + '\t'
        let out_2 = (avg_biomass_t.toExponential(2)) + '\t'
        let out_3 = (avg_biomass_c.toExponential(2)) + '\n'
        let out_0 = this.time + '\t'
        fs.appendFileSync(resource_count, out_0 + out_1 + out_2 + out_3)

        let out_4 = (avg_Vmax_b.toExponential(2)) + '\t'
        let out_5 = (avg_Vmax_t.toExponential(2)) + '\n'
        fs.appendFileSync(Vmax_file, out_0)
        fs.appendFileSync(Vmax_file, out_4)
        fs.appendFileSync(Vmax_file, out_5)
    }

    if (this.time % mutation_data_interval === 0) {
        let out_1 = this.time + '\t'
        let out_2 = b_to_c + '\t'
        let out_3 = c_to_b + '\t'
        let out_4 = t_to_c + '\t'
        let out_5 = c_to_t + '\n'
        fs.appendFileSync(mutation_count_vmax, out_1)
        fs.appendFileSync(mutation_count_vmax, out_2)
        fs.appendFileSync(mutation_count_vmax, out_3)
        fs.appendFileSync(mutation_count_vmax, out_4)
        fs.appendFileSync(mutation_count_vmax, out_5)
        
        b_to_c = 0;
        c_to_b = 0;
        t_to_c = 0;
        c_to_t = 0;
    }

    if (sim.time % 10000 === 0) {       // logs lineages every 10000 time steps
        for (let i = 0; i < this.nc; i++) {
            for (let j = 0; j < this.nr; j++) {
                let cell = this.grid[i][j];
                if (cell.type === 'b') {
                    let uID = cell.uID;
                    uid_distribution_b[uID] = (uid_distribution_b[uID] || 0) + 1;
                } else if (cell.type === 't') {
                    let uID = cell.uID;
                    uid_distribution_t[uID] = (uid_distribution_t[uID] || 0) + 1;
                }
            }
        }
    
        // Categories for pioneers and harvesters, including special categories for p_c and h_c
        const categories_b = [
            "p_c",
            "0 to 1",
            "1 to 2",
            "2 to 3",
            "3 to 4",
            "4 to 5",
            "5 to 6",
            "6 to 7",
            "7 to 8",
            "8 to 9",
            "9 to 10",       // For values greater than or equal to 9
            "10 to 11",
            "11 to 12",
            "12 to 13",
            "13 to 14",
            "14 to 15",
            "15 to 16",
            "16 to 17",
            "17 to 18",
            "18 to 19",
            "19 to 20",
            "20 to 21",
            "21 to 22",
            "22 to 23",
            ">23"
        ];
    
        const categories_t = [
            "h_c", 
            "0 to 1",
            "1 to 2",
            "2 to 3",
            "3 to 4",
            "4 to 5",
            "5 to 6",
            "6 to 7",
            "7 to 8",
            "8 to 9",
            "9 to 10",
            "10 to 11",
            "11 to 12",
            "12 to 13",
            "13 to 14",
            "14 to 15",
            "15 to 16",
            "16 to 17",
            "17 to 18",
            "18 to 19",
            "19 to 20",
            "20 to 21",
            "21 to 22",
            "22 to 23",
            ">23"
        ];
    
        // Initialize category counts for pioneers and harvesters
        let categoryCountsB = {};
        let categoryCountsT = {};
    
        categories_b.forEach(category => categoryCountsB[category] = 0);
        categories_t.forEach(category => categoryCountsT[category] = 0);
    
        // Categorize uIDs for pioneers
        for (let uID in uid_distribution_b) {
            let category = categories_b.find(c => uID === c || uID.includes(c));
            if (category) {
                categoryCountsB[category] += uid_distribution_b[uID];
            }
        }
    
        // Categorize uIDs for harvesters
        for (let uID in uid_distribution_t) {
            let category = categories_t.find(c => uID === c || uID.includes(c));
            if (category) {
                categoryCountsT[category] += uid_distribution_t[uID];
            }
        }
    
        // Writing CSV content or logging for pioneers and harvesters
        let csvContentB = 'Category,Count\n';
        let csvContentT = 'Category,Count\n';
    
        categories_b.forEach(category => {
            csvContentB += `${category},${categoryCountsB[category]}\n`;
        });
    
        categories_t.forEach(category => {
            csvContentT += `${category},${categoryCountsT[category]}\n`;
        });
    
        // Append to respective files
        fs.appendFileSync(lineage_dist_b, csvContentB);
        fs.appendFileSync(lineage_dist_t, csvContentT);
    }
    

    if (sim.time % dist_interval === 0) {       // logs trait value distribution 
        for (let i = 0; i < this.nc; i++) {
            for (let j = 0; j < this.nr; j++) {
                if (this.grid[i][j].type === 'b') {
                    vmax_value_b = (this.grid[i][j].Vmax_polymer) * (Math.pow(PG_size, 2) / dT)
                    vmax_dist_b.push(vmax_value_b)
                }
                else if (this.grid[i][j].type === 't') {
                    vmax_value_t = (this.grid[i][j].Vmax_polymer)/dT
                    //console.log(vmax_value_t)
                    vmax_dist_t.push(vmax_value_t)
                }    
            }
        }
       
        // Define your fixed bin ranges
        const binRanges = [
            "0",       // For values equal to 0 (this is noise that should be removed)
            "0 - 1",
            "1 - 2",
            "2 - 3",
            "3 - 4",
            "4 - 5",
            "5 - 6",
            "6 - 7",
            "7 - 8",
            "8 - 9",
            "9 - 10",       // For values greater than or equal to 9
            "10 - 11",
            "11 - 12",
            "12 - 13",
            "13 - 14",
            "14 - 15",
            "15 - 16",
            "16 - 17",
            "17 - 18",
            "18 - 19",
            "19 - 20",
            "20 - 21",
            "21 - 22",
            "22 - 23",
            ">23",
        ];
        

        // Initialize bin counts
        let binCountsB = Array(binRanges.length).fill(0);
        let binCountsT = Array(binRanges.length).fill(0);
        
        // Calculate bin indices and update bin counts based on your fixed bin ranges
        vmax_dist_b.forEach((vmax_b) => {
            if (!isNaN(vmax_b)) { // Check if vmax is a valid number
                let binIndex = findBinIndex(vmax_b); // Implement a function to find the appropriate bin index
                binCountsB[binIndex]++;
            }
        });

        vmax_dist_t.forEach((vmax_t) => {
            if (!isNaN(vmax_t)) {
                let binIndex = findBinIndex(vmax_t);
                binCountsT[binIndex]++;
            }
        });

        function findBinIndex(value) {
            for (let i = 0; i < binRanges.length; i++) {
                const range = binRanges[i];
                const parts = range.split('-');
                if (parts.length === 1) {
                    // Single value range, e.g., "0" or ">23"
                    const binValue = parseFloat(parts[0]);
                    if (value === binValue || (parts[0] === '>' && value > binValue)) {
                        return i;
                    }
                }     
                else if (parts.length === 2) {
                    // Range with two values, e.g., "0 - 1"
                    const minBinValue = parseFloat(parts[0]);
                    const maxBinValue = parseFloat(parts[1]);
                    if (value > minBinValue && value <= maxBinValue) {
                        return i;
                    }
                }
            }
            return binRanges.length - 1;
        }    
        
        // Modify the distribution data based on the fixed bin ranges
        const distributionDataB = binRanges.map((range, index) => ({
            Bin: range.replace(/\s-\s/, ' to '), // Adjust bin label formatting
            Count: binCountsB[index],
        }));

        const distributionDataT = binRanges.map((range, index) => ({
            Bin: range.replace(/\s-\s/, ' to '),
            Count: binCountsT[index],
        }));

        // Process the distribution data
        distributionDataB.forEach((data) => {
            appendDistributionData(data, vmax_dist_b);
        });

        distributionDataT.forEach((data) => {
            appendDistributionData(data, vmax_dist_t);
        });
    }

    let total_polymer_digested = 0
    if (poly_digestion_sum.length > 0) {
        total_polymer_digested = poly_digestion_sum.reduce(function(a, b) {
        return a + b;
        });
        master_poly_digestion_sum.push(total_polymer_digested);
        //console.log(total_polymer_digested)
        poly_digestion_sum.length = 0
    }
}

sim.resources.update = function() {
    this.time++;
    if (mix === false) {
        this.diffuseStates('monomer_count', diffusion_rate * dT)
    }

    let sum_polymer = 0
    let sum_monomer = 0
    let sum_monomer_b = 0
    let sum_b = 0
    let avg_b = 0
	let sum_monomer_t = 0
	let sum_t = 0
	let avg_t = 0
    let sum_monomer_c = 0
    let sum_c = 0
    let avg_c = 0
                    
    for (let i = 0; i < this.nc; i++) {
        for (let j = 0; j < this.nr; j++) {
            sum_polymer += this.grid[i][j].polymer_count    // counts the total polymer count in the grid
            sum_monomer += this.grid[i][j].monomer_count    // counts the total monomer count in the grid
            if (sim.colonies.grid[i][j].type === 'b') {
                sum_b++;
                sum_monomer_b += this.grid[i][j].monomer_count
            }
			else if (sim.colonies.grid[i][j].type === 't') {
                sum_t++;
                sum_monomer_t += this.grid[i][j].monomer_count
            }
            else if (sim.colonies.grid[i][j].type === 'c') {
                sum_c++;
                sum_monomer_c += this.grid[i][j].monomer_count
            }
        }                    
    }

    avg_b = sum_monomer_b / sum_b
    avg_t = sum_monomer_t / sum_t
    avg_c = sum_monomer_c / sum_c

    if (this.time % data_interval === 0) {
        let out_1 = (avg_b.toExponential(2)) + '\t'
		let out_2 = (avg_t.toExponential(2)) + '\t'
        let out_3 = (avg_c.toExponential(2)) + '\n'
        let out_0 = this.time + '\t'
        let out_m = (sum_monomer.toExponential(2)) + '\n'
        fs.appendFileSync(bias_file, out_0 + out_1 + out_2 + out_3)
        fs.appendFileSync(monomer_count, out_0 + out_m)
    }
    
    if (sum_polymer < resource_threshold) {             // check if particle is exhausted for a serial transfer
        let master_poly_sum = 0
        if (master_poly_digestion_sum.length > 0) {
            master_poly_sum = master_poly_digestion_sum.reduce(function(a, b){return a + b});
            //console.log(master_poly_sum)
            master_poly_digestion_sum.length = 0;
        }    
        if (b_neck === true) {
            sim.colonies.bottleneck(0.9)            // serial transfer triggered
        }
        bottleneck_count++;   
        sim.resources.refresh() 
        let time_between_bottlenecks = this.time - bottleneck_time;
        bottleneck_time = this.time;
        let out_1 = bottleneck_time + '\t'
        let out_2 = time_between_bottlenecks + '\t'
        let out_3 = master_poly_sum.toExponential(2) + '\t'
        let out_4 = (out_3/out_2).toExponential(2) + '\n'
        fs.appendFileSync(bottleneck_time_file, out_1 + out_2 + out_3 + out_4)
        
        sim.colonies.perfectMix()               // mix the remaining cells so that they are randomly placed on the fresh particle
    }

    // for generating png images or txt heatmaps

    if (run_num === 1 && this.time < 15000) {           // i opted for a few images up to a certain time but this can be changed as needed
        if (bottleneck_count / image_duration === 5) {  // for example, here i want to trigger image capture from the 5th serial transfer
            grid_capture = true;
        }
           
        if (grid_capture === true) {            // i generated images from when a serial transfered occurred all the way up to another serial transfer
            if (bottleneck_count === image_duration + 2) {      // the prorgam will stop taking images upon triggering the 7th serial transfer
                grid_capture = false;
                bottleneck_count = 0;                   // reset the serial transfer counter, and start taking images again at the 5th one
                //console.log("done")
            };
    
            if (take_png === true) {
                for (let i = 0; i < sim.colonies.nc; i++) {
                    for (let j = 0; j < sim.colonies.nr; j++) {
                        let cell = sim.colonies.grid[i][j]
                        color = '#000000'
                        if (cell.type === 'b') {
                            color = '#9F2B68';  // amaranth
                        }
                        else if (cell.type === 't') {
                            color = '#FFA700'; // cyan
                        }
                        else if (cell.type === 'c') {
                            color = '#00FFFF'; // cyan
                        }
                        ctx.fillStyle = color   
                        ctx.fillRect(i*cellSize, j*cellSize, cellSize, cellSize)
                    }
                }
                buf = canvas.toBuffer('image/png');
                writeFileSync(join(movie, `grid_${this.time}.png`), b)
                ctx.clearRect(0, 0, canvasWidth, canvasHeight);    // clear the canvas after every image genration to save RAM
                //sim.resources.heatmap_polymer();
                buf = null;  // to mark the buffer for garbage collection
            }
            else if (take_txt === true) {
                //console.log(take_txt+"inside")
                sim.resources.heatmap_monomer();
                sim.colonies.heatmap_type_biomass();
            } 
        }
    };

    if (sum_b === 0 || sum_t === 0 || sum_c === 0) {
        if (!extinction_event) {
            extinction_event = true;
            extinction_time = this.time
        }
    };
    
    if (extinction_event && this.time > extinction_time + 10000) {      // choose when to end the simulation, upon exclusion of a cell type, here
        process.exit()
    }
}

//sim.addButton("Play/Pause", function () { sim.toggle_play() })

sim.initialise()

sim.start()
