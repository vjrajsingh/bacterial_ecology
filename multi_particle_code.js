// multi-particle code - broadcasters and tetherers only

/*  the following is an attempt to create a multi-particle simulation using cacatoo.
    to my knowledge, and from what my co-supervisor has told me, this has not been done before,
    so let us see what we can do! */

// node js requirements

let Simulation = require('cacatoo') // Loads the Simulation class from local package like below
let yargs = require('yargs')
let fs = require('fs')              // for creating, editing and appending files
let MersenneTwister = require('mersenne-twister') // for seeding
const { createCanvas } = require('@napi-rs/canvas')
const { join } = require('path')
const { writeFileSync } = require('fs')

//const readline = require('readline'); // reads data input from the cosnole

let numGrids = 30 // how many grids/particles? default value 30

//console.log(gridIndices)

function shuffle(array) {
    for (let q = array.length - 1; q > 0; q--) {
        let r = Math.floor(Math.random() * (q + 1));
        [array[q], array[r]] = [array[r], array[q]];
    }
    return array
}

let death = 0.1
let birth = 0.5
let migration                        // since this parameter is being scanned, there is no default value
let destruction                      // same for this   
let PG_size = 3
let maxsteps
let num_col = 100 // per grid
let num_row = 100 // per grid
let starting_polymer_value          // and for this
let data_interval = 100
let div_threshold = 2
//let resource_threshold = 1000
let starting_resource_value = 1.0
let broadcaster_fraction = 0.1
let tetherer_fraction = 0.1
let cheater_fraction = 0.1
let extinction_time = 0.0
let extinction_event = false
let grid_capture = false
let grid_index_for_images;
let image_interval = 300;

let dT = 0.1 // [0, 1]
//let b = 0.3 // percentage of population 'lost' after grid exhaustion
let mt = new MersenneTwister(); // will create a new MersenneTwister instance

// timestamps of interest //

let timestamps_of_interest = [100, 100000, 950000]
let time_at_grid_capture = 0;

// all time-related variables have to be divided by dT

let migration_time = 100; 
let destroyed_refresh = 1;
let diff_interval = 1; 

//bacterial properties

let MM_polymers = 0.5
let MM_monomers = 0.2
let aly_b = 9
let aly_t = 0.5
let oal_all = 0.5
let aly_cost = 0.01
let oal_cost = 0.01
let migration_cost = 0.01

let n = mt.random_int() // 32-bit range. has highest entropy

let run_num

// block for parsing arguments 

const { hideBin } = require('yargs/helpers')
let argv = yargs(hideBin(process.argv)).argv

if (typeof argv.numGrids !== "undefined") {
	numGrids = argv.numGrids
	console.log("Number of grids in the system -- user specified", numGrids)
} 
else {
	console.log("Number of grids in the system\t", numGrids)
}

if (typeof argv.death !== "undefined") {
	death = argv.death
	console.log("Random Death -- user specified", death)
} 
else {
	console.log("Random Death\t", death)
}

if (typeof argv.migration !== "undefined") {
	migration = argv.migration
	console.log("Migration rate -- user specified", migration)
} 
else {
	console.log("Migration rate\t", migration)
}

if (typeof argv.destruction !== "undefined") {
	destruction = argv.destruction
	console.log("Random grid destruction rate -- user specified", destruction)
} 
else {
	console.log("Random grid destruction rate\t", destruction)
}

if (typeof argv.aly_cost !== "undefined") {
	aly_cost = argv.aly_cost
	console.log("Loss due to aly secretion -- user specified", aly_cost)
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

if (typeof argv.migration_cost !== "undefined") {
	migration_cost = argv.migration_cost
	console.log("Loss due to migration", migration_cost)
} 
else {
	console.log("Loss due to migration\t", migration_cost)
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
	console.log("Aly enzyme value for tetherer\t", aly_t)
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
    console.log("Km value for polymer digestion -- user specified", MM_polymers)
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

/*if (typeof argv.resource_threshold !== "undefined") {
	resource_threshold = argv.resource_threshold
	console.log("Critical resource value after which resources refresh -- user specified", resource_threshold)
} 
else {
	console.log("Critical resource value after which resources refresh\t", resource_threshold)
}*/

if (typeof argv.starting_polymer_value !== "undefined") {
	starting_polymer_value = argv.starting_polymer_value
	console.log("Initial polymer count per grid point -- user specified", starting_polymer_value)
} 
else {
	console.log("Initial polymer count per grid point\t", starting_polymer_value)
}

if (typeof argv.run_num !== "undefined") {
    run_num = argv.run_num
    console.log("Run number -- user specified", run_num)
}
else {
    console.log("Run number\t", run_num)
}

if (typeof argv.migration_time !== "undefined") {
    migration_time = argv.migration_time
    console.log("Time interval between migrations -- user specified", migration_time)
}
else {
    console.log("Time interval between migrations\t", migration_time)
}

if (typeof argv.destroyed_refresh !== "undefined") {
    destroyed_refresh = argv.destroyed_refresh
    console.log("Randomly destroyed grid replacement rate -- user specified", destroyed_refresh)
}
else {
    console.log("Randomly destroyed grid replacement speed\t", destroyed_refresh)
}

if (typeof argv.diff_interval !== "undefined") {
    diff_interval = argv.diff_interval
    console.log("Time interval between margolus diffusion or well-mixing -- user specified", diff_interval)
}
else {
    console.log("Time interval between margolus diffusion or well-mixing\t", diff_interval)
}

// all of my boolean variables I want passed as comamnd line arguments (without any JS bullshit of type coercion) are handled here

// one can choose to take a snapshot of one randomly chosen grid during runtime using the node API canvas module. there is a buffer issue here as the canvas image gets pasted to a buffer object and then converted to a png file from the binary data. this is causing memory issues and if someone can fix that, all the more better!
// or they can convert the grid to a txt file and later convert that to an image

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

// Now you can access these options as properties of the argv object
let mar_diffusion = argv.mar_diffusion;
let mix = argv.mix;
let take_txt = argv.take_txt;
let take_png = argv.take_png;

console.log('Margolus Diffusion?:', mar_diffusion);
console.log('Perfectly Mixed?:', mix);
console.log('Making text files of the grid?:', take_txt);
console.log('Making png files of the grid in runtime?:', take_png);

let gridIndices = Array.from({ length: numGrids }, (_, index) => index); // convert the grid indices into an array of grid indices

let home = 'output/poly_'+starting_polymer_value+'mig_'+migration+'dest_'+destruction
if (!fs.existsSync(home)){
    fs.mkdirSync(home);
}

let dir = home+'/data_out_'+run_num
if (!fs.existsSync(dir)){
    fs.mkdirSync(dir);
}

let pop = dir+"/population_data"
if (!fs.existsSync(pop)) {
    fs.mkdirSync(pop)
}

let migrants = dir+"/migrants_data"
if (!fs.existsSync(migrants)) {
    fs.mkdirSync(migrants)
}

let pop_output_files = [];                 // to store population data files
//let internal_resource_output_files = [];   // to store biomass data files
let migrants_output_files = [];            // to store relevant migrant information data files
//let migration_counters_output_files = [];  // to store the 'count' of migrations occurring per grid
//let avg_monomer_output_files = []


// track the destruction of grids?
// track the destruction of grids?
let destroyed_grids_file = dir+"/destroyed_grids.dat"    // stores the grid number and its destruction time stamp
if (fs.existsSync(destroyed_grids_file)) {
    fs.unlinkSync(destroyed_grids_file)
}

for (let i = 0; i < numGrids; i++) {
    let main_file = pop +`/output_${i}.dat`
    if (fs.existsSync(main_file)) {
        fs.unlinkSync(main_file)
    }
    pop_output_files.push(main_file)

    let migrants_file = migrants +`/migrations_details_${i}.dat`
    if (fs.existsSync(migrants_file)) {
        fs.unlinkSync(migrants_file)
    }
    migrants_output_files.push(migrants_file) 
}

let master_main_file = dir+"/system_population_size.dat"   // stores global population size
if (fs.existsSync(master_main_file)) {
    fs.unlinkSync(master_main_file)
}

let master_biomass_file = dir+"/system_internal_resources.dat"  // stores global biomass info
if (fs.existsSync(master_biomass_file)) {
    fs.unlinkSync(master_biomass_file)
}

let destruction_count = dir+"/total_disturbances.dat"    // outputs an array which is a counter for each grids destruction
if (fs.existsSync(destruction_count)) {
    fs.unlinkSync(destruction_count)
}

let movie;
let heatmaps;

if (take_png === true && run_num === 1) {
    movie = dir+"/images"
    if (!fs.existsSync(movie)) {
        fs.mkdirSync(movie)
    }
}    

if (take_txt === true && run_num === 1) {
    heatmaps = dir+"/heatmaps"
    if (!fs.existsSync(heatmaps)) {
        fs.mkdirSync(heatmaps)
    }
}

let run_info = dir+"/Run_Info.txt"
if(fs.existsSync(run_info)){
//	console.log("in exists block\t\t", filn)
	fs.unlinkSync(run_info)
}

let final_result = `output/results/result_run_${run_num}_poly_${starting_polymer_value}_mig_${migration}_dest_${destruction}.json`
/*if (fs.existsSync(final_result)) {
    fs.unlinkSync(final_result)
}*/

let destroyed_refresh_final = destroyed_refresh / dT
let diff_interval_final = diff_interval / dT

let aly_broadcaster = aly_b * dT
let aly_tetherer = aly_t * dT
let oal = oal_all * dT

// all 'rates'
let rand_death = death * dT
let f = birth * dT
let migration_rate = migration * dT
let destruction_rate = destruction * dT

//console.log( migration, migration_rate)

fs.appendFileSync(run_info, "dT =\t"+dT+"\t"+"All values are ignoring dT")
fs.appendFileSync(run_info, "Run number =\t"+run_num+"\n"+"Public goods neighbourhood size =\t"+PG_size+"\n"+"Random death probability =\t"+death+"\n"+"Fecundity rate =\t"+birth+"\n"+"Migration rate =\t"+migration+"\n"+"Random grid destruction rate =\t"+destruction+"\n")
fs.appendFileSync(run_info, "Maxtime =\t"+maxsteps+"\n"+"Grid size =\t"+num_col+" by "+num_row+"\n")
fs.appendFileSync(run_info, "Cost of digesting polymers =\t"+aly_cost+"\n"+"Cost of digesting monomers in broadcasters =\t"+oal_cost+"\n"+"Cost of migration =\t"+migration_cost+"\n"+"Division threshold =\t"+div_threshold+"\n")
fs.appendFileSync(run_info, "Enzyme value for broadcasters =\t"+aly_b+"\n"+"Enzyme value for tetherers =\t"+aly_t+"\n"+"Enzyme value for monomer digestion (all cells) =\t"+oal_all+"\n")
fs.appendFileSync(run_info, "Michaelis-Menten constant (polymer) =\t"+MM_polymers+"\n"+"Michaelis-Menten constant (monomer) =\t"+MM_monomers+"\n")
fs.appendFileSync(run_info, "Perfectly mixed? =\t"+mix+"\n"+"Margolus Diffusion? =\t"+mar_diffusion+"\n"/*+"Public Goods Operon? =\t"+PG_switch+"\n"*/)
fs.appendFileSync(run_info, "Write out interval =\t"+data_interval+"\n"+"Image capture interval =\t"+image_interval+"\n"/*+"Resource refresh interval =\t"+resource_threshold+"\n"*/+"Polymer units per grid point =\t"+starting_polymer_value+"\n")
fs.appendFileSync(run_info, "Time interval between migrations =\t"+migration_time+"\n"+"Grid recovery rate after destruction =\t"+destroyed_refresh+"\n"+"Time interval between diffusions or mixing =\t"+diff_interval+"\n")
fs.appendFileSync(run_info, "Ratio of broadcaster to tetherer to cheater =\t"+broadcaster_fraction+":"+tetherer_fraction+":"+cheater_fraction+"\t")
fs.appendFileSync(run_info, "Total number of grids in the system =\t"+numGrids+"\t"+"Taking png? =\t"+take_png+"\n"+"Taking txt? =\t"+take_txt+"\n")

let config = {
    maxtime: maxsteps,
    ncol: num_col,
    nrow: num_row,
    seed: n,
    wrap: [true, true],
    scale: 3,
    statecolours: {'type':{'b': 'red', 't': 'yellow', 'c': 'green'}}
}

fs.appendFileSync(run_info, "seed =\t"+config.seed+"\n")

sim = new Simulation(config)

// create array to hold references to colony layers (grid models).
let colonyGridModels = [];
let resourcesGridModels = [];
// array to hold the destroyed grids
let destroyedGrids = [];
// make arrays for your grid's occupants
let broadcaster_populations = new Array(numGrids).fill(0);
let tetherer_populations = new Array(numGrids).fill(0);
let cheater_populations = new Array(numGrids).fill(0);
let broadcaster_internal_resources = new Array(numGrids).fill(0);
let tetherer_internal_resources = new Array(numGrids).fill(0);
let cheater_internal_resources = new Array(numGrids).fill(0);
let destruction_counters = new Array(numGrids).fill(0);

// system-wide letiables
let total_broadcaster_system = 0;
let total_tetherer_system = 0;
let total_cheater_system = 0;
let total_internal_resources_system_p = 0;
let total_internal_resources_system_h = 0;
let total_internal_resources_system_c = 0;

// create a list of properties for your microbes
let microbes = [{type: 'b', Vmax_polymer: aly_broadcaster/(Math.pow(PG_size, 2)), Km_polymer: MM_polymers, Vmax_monomer: oal, Km_monomer: MM_monomers, internal_resources: starting_resource_value},
                {type: 't', Vmax_polymer: aly_tetherer, Km_polymer: MM_polymers, Vmax_monomer: oal, Km_monomer: MM_monomers, internal_resources: starting_resource_value},
                {type: 'c', Vmax_polymer: 0.0, Km_polymer: 0.0, Vmax_monomer: oal, Km_monomer: MM_monomers, internal_resources: starting_resource_value}]

// loop through each grid layer.
for (let i = 0; i < numGrids; i++) {
    const colonyModelName = `colonies_${i}`;
    
    // create a grid model with a unique name (my e.g., "colonies_1", "colonies_2", etc.).
    sim.makeGridmodel(colonyModelName);
    
    /* how many grids should be colonised beforehand? you can use this block to test that
    if (i < 5) {
        sim.populateGrid(colonyModelName, microbes, [broadcaster_fraction, tetherer_fraction, cheater_fraction]); // initialised grids
    } else {
        sim.populateGrid(colonyModelName, microbes, [0.0, 0.0, 0.0]); // the others are not initialised
    }*/

    // or choose to initialise the whole system with your microbes
    sim.populateGrid(colonyModelName, microbes, [0.0, 0.0, 0.0]);         // populateSpot() returns undefined for empty grid points and to avoid that, we initialise every grid point to 0 beforehand
    sim.populateSpot(colonyModelName, microbes, [broadcaster_fraction, tetherer_fraction, cheater_fraction], 2000, config.ncol/2, config.nrow/2);

    // push the grid model into the array of grid models.
    colonyGridModels.push(sim[colonyModelName]);

    // create displays
    //sim.createDisplay(colonyModelName, 'type', "");
    
    //console.log(colonyGridModels)
    
    // same for the resources grid models
    
    const resourcesModelName = `resources_${i}`;
    
    sim.makeGridmodel(resourcesModelName);
    
    sim.initialGrid(resourcesModelName, 'polymer_count', starting_polymer_value, 1.0)
    sim.initialGrid(resourcesModelName, 'monomer_count', 0.0, 1.0)
    
    resourcesGridModels.push(sim[resourcesModelName]);
    
    /*
    sim.createDisplay_continuous({
    model: resourcesModelName, // Use a unique model name for each layer
    property: 'polymer_count',
    label: `Polymer count Layer ${i}`,
    minval: 0,
    maxval: 3,
    fill: 'viridis'
    });
    
    
    sim.createDisplay_continuous({
    model: resourcesModelName, // Use a unique model name for each layer
    property: 'monomer_count',
    label: `Monomer count Layer ${i}`,
    minval: 0,
    maxval: 3,
    fill: 'inferno'
    });
    */
   // console.log(resourcesModelName) 
}

// in multi-particle system, one has to define the next state and update logics separately in a function. this function contains several nested functions that will be applied to each
// subsequent grid model

// each grid must have a unique model number, which will be used to track how many times it has been destroyed. for example, if '3 - 5', it means grid 3 has been destroyed 5 times.
for (let x = 0; x < resourcesGridModels.length; x++) {
    let resources = resourcesGridModels[x];
    //let colonies = colonyGridModels[x];
    resources.gridNumber = x;
    resources.refreshCounter = 0;
    resources.uID = `${resources.gridNumber}-${resources.refreshCounter}`;
    //colonies.imageCaptureState = false
    //console.log(resources.uID)
}

// using node's canvas module to capture grid images (optional)
const canvasWidth = num_row * config.scale;
const canvasHeight = num_col * config.scale;
const cellSize = config.scale; 
let canvas = createCanvas(canvasWidth, canvasHeight);
let ctx = canvas.getContext('2d');
let b;
let color;

function defineNextStates() {
    for (let x = 0; x < gridIndices.length; x++) {   // the number of grids in the system is the same as the number of elements in the gridIndices array
        let colonies = colonyGridModels[x];
        let resources = resourcesGridModels[x];
        //console.log(colonies)
        //console.log(resources)
        //console.log(x)

        colonies.nextState = function(i, j) {               
            //console.log(grid_index)
            let randomneigh = colonies.randomMoore8(this, i, j)
            let this_gp = colonies.grid[i][j]
            let resources_gp = resources.grid[i][j]

        // chance to die either randomly or due to lack of resources   

            if (this_gp.internal_resources < 0.0) {     // death due to lack of internal resources
                this_gp.type = 0.0
                this_gp.Km_polymer = 0.0
                this_gp.Vmax_polymer = 0.0
                this_gp.Km_monomer = 0.0
                this_gp.Vmax_monomer = 0.0
                this_gp.internal_resources = 0.0
            }

            if (this.rng.genrand_real1() < rand_death) { // random death every time step
                this_gp.type = 0.0
                this_gp.Km_polymer = 0.0
                this_gp.Vmax_polymer = 0.0
                this_gp.Km_monomer = 0.0
                this_gp.Vmax_monomer = 0.0
                this_gp.internal_resources = 0.0
            }

            // reproduction loop

            if (!this_gp.type && randomneigh.internal_resources > div_threshold) { // reproduce if you have enough resources. all mutation logic needs to be applied here
                if (this.rng.genrand_real1() < f) {
                    this_gp.type = randomneigh.type
                    this_gp.Km_polymer = randomneigh.Km_polymer
                    this_gp.Vmax_polymer = randomneigh.Vmax_polymer
                    this_gp.Km_monomer = randomneigh.Km_monomer
                    this_gp.Vmax_monomer = randomneigh.Vmax_monomer
                    randomneigh.internal_resources = this_gp.internal_resources =  randomneigh.internal_resources / 2
                }
            }

            // monomer digestion loop starts here    

            if (resources_gp.monomer_count > 0.0) {
                switch (this_gp.type) {
                    case 'c':
                        let mono_digestion_c = (this_gp.Vmax_monomer * resources_gp.monomer_count) / (this_gp.Km_monomer + resources_gp.monomer_count)
                        if (resources_gp.monomer_count < mono_digestion_c) {resources_gp.monomer_count = mono_digestion_c}
                        resources_gp.monomer_count -= mono_digestion_c
                        this_gp.internal_resources += mono_digestion_c
                        this_gp.internal_resources -= oal_cost * this_gp.Vmax_monomer
                    break;

                    case 't':
                        let mono_digestion_h = (this_gp.Vmax_monomer * resources_gp.monomer_count) / (this_gp.Km_monomer + resources_gp.monomer_count)
                        if (resources_gp.monomer_count < mono_digestion_h) {resources_gp.monomer_count = mono_digestion_h}
                        resources_gp.monomer_count -= mono_digestion_h
                        this_gp.internal_resources += mono_digestion_h
                        this_gp.internal_resources -= oal_cost * this_gp.Vmax_monomer
                    break;

                    case 'b':
                        let mono_digestion_p = (this_gp.Vmax_monomer * resources_gp.monomer_count) / (this_gp.Km_monomer + resources_gp.monomer_count)
                        if (resources_gp.monomer_count < mono_digestion_p) {resources_gp.monomer_count = mono_digestion_p}
                        resources_gp.monomer_count -= mono_digestion_p
                        this_gp.internal_resources += mono_digestion_p
                        this_gp.internal_resources -= oal_cost * this_gp.Vmax_monomer
                    break;
                }
            }


        // polymer digestion loop starts here

            if (resources_gp.polymer_count > 0.0 && this_gp.type === 't') {
                let poly_digestion_h = (this_gp.Vmax_polymer * resources_gp.polymer_count) / (this_gp.Km_polymer + resources_gp.polymer_count)
                if (resources_gp.polymer_count < poly_digestion_h) {resources_gp.polymer_count = poly_digestion_h}
                resources_gp.polymer_count -= poly_digestion_h
                resources_gp.monomer_count += poly_digestion_h
                this_gp.internal_resources -= aly_cost * this_gp.Vmax_polymer
            }

            else if (resources_gp.polymer_count > 0.0 && this_gp.type === 'b') {
                let poly_digestion_p = []
                switch (PG_size) {
                    case 3:                                         // public goods neighbourhood of 3
                        for (let n = 0; n <= 8; n++) {
                            poly_digestion_p[n] = (this_gp.Vmax_polymer * this.getNeighbour(resources, i, j, n).polymer_count) / (this_gp.Km_polymer + this.getNeighbour(resources, i, j, n).polymer_count)
                            if (this.getNeighbour(resources, i, j, n).polymer_count < poly_digestion_p[n]) {this.getNeighbour(resources, i, j, n).polymer_count = poly_digestion_p[n]}
                            this.getNeighbour(resources, i, j, n).polymer_count -= poly_digestion_p[n]    // access the neighbourhood of the gridpoint and update the properties accrodingly. 
                            this.getNeighbour(resources, i, j, n).monomer_count += poly_digestion_p[n]
                        }
                        //console.log(poly_digestion_p)
                        this_gp.internal_resources -= (aly_cost * this_gp.Vmax_polymer * Math.pow(PG_size, 2)) // cost of public goods production is constant cost times the enzyme units secreted.
                    break; 
                }
            }
        }

        resources.nextState = function(i, j) {}
    }      
}      

// here is a function of nested functions where we define cacatoo's update functions to be applied to the entire grid per time step
function defineUpdates() {
    for (let x = 0; x < colonyGridModels.length; x++) {
        let colonies = colonyGridModels[x];
        let resources = resourcesGridModels[x];
        //console.log(colonies)
        //console.log(resources)
        //console.log(0)
        let destruction_counter = 0;
        let migration_counter = 0;
        let invasion_counter = 0;

        function refreshDestroyedGrids() {
            // Get the current simulation time
            let refresh_time = sim.time;

            for (let i = 0; i < destroyedGrids.length; i++) {
                let destroyedGrid = destroyedGrids[i];
                let destroyedResource = destroyedGrid.resource;
                let destroy_time = destroyedGrid.destroy_time;

                // Check if enough time steps have passed since destruction
                    
                if (refresh_time === destroy_time + destroyed_refresh_final) {
                    // Reset resource grid
                    //console.log("done")
                    sim.initialGrid(destroyedResource, 'polymer_count', starting_polymer_value, 1.0);
                    sim.initialGrid(destroyedResource, 'monomer_count', 0.0, 1.0);
                    //console.log(resources)

                    destroyedResource.refreshCounter++;
                    destroyedResource.uID = `${destroyedResource.gridNumber}-${destroyedResource.refreshCounter}`;
                    //console.log(destroyedResource.uID)
                    //console.log(resources)
                    //console.log(destroyedResource, destroyedResource.destruction_counter)
                    
                    // Remove this grid from the list since it has been refreshed
                    destroyedGrids.splice(i, 1);
                    i--; // Decrement i to account for the removed element
                }
            }
        }
        
        // most important function (for me since this was a pain to implement) since cacatoo's original step() function iterates over all the gridmodels in sequential order,
        // it was necessary to overwrite the step function. since there is no parallel execution, we will update the grids in a random order. the array of grids (gridIndices)
        // is shuffled using the Fisher-Yates algorithm and then the step function is modified as follows:

        sim.step = function() {
            shuffle(gridIndices)           // shuffle the gridIndices array
            //console.log(gridIndices)

            for (let i of gridIndices) {    // then implement the update sequence based on the shuffled array
                colonyGridModels[i].update();
                resourcesGridModels[i].update();
            }
            sim.time++;                    // increment time so that timesteps increase ONLY after ALL the grids have been updated
            total_broadcaster_system = broadcaster_populations.reduce(function(a, b){return a + b});
            total_tetherer_system = tetherer_populations.reduce(function(a, b){return a + b});
            total_cheater_system = cheater_populations.reduce(function(a, b){return a + b});
            total_internal_resources_system_p = broadcaster_internal_resources.reduce(function(a, b){return a + b});
            total_internal_resources_system_h = tetherer_internal_resources.reduce(function(a, b){return a + b});
            total_internal_resources_system_c = cheater_internal_resources.reduce(function(a, b){return a + b});

            //console.log("ugh")

            if (sim.time % data_interval === 0) {
                let out_1 = sim.time + '\t'
                let out_2 = total_broadcaster_system + '\t'
                let out_3 = total_tetherer_system + '\t'
                let out_4 = total_cheater_system + '\n'
                let out_5 = (total_internal_resources_system_p/total_broadcaster_system).toExponential(2) + '\t'
                let out_6 = (total_internal_resources_system_h/total_tetherer_system).toExponential(2) + '\t'
                let out_7 = (total_internal_resources_system_c/total_cheater_system).toExponential(2) + '\n'

                fs.appendFileSync(master_main_file, out_1)
                fs.appendFileSync(master_main_file, out_2)
                fs.appendFileSync(master_main_file, out_3)
                fs.appendFileSync(master_main_file, out_4)
                fs.appendFileSync(master_biomass_file, out_1)
                fs.appendFileSync(master_biomass_file, out_5)
                fs.appendFileSync(master_biomass_file, out_6)
                fs.appendFileSync(master_biomass_file, out_7)
            }

            /* reworking this section as there are three population types now
               the program will stop when any one population type goes extinct
               if broadcasters are the first to die, print 1, and then if everyone dies, print 0 as well, else print 1
               if tetherers are the first to die, print 2, and then if everyone dies, print 0 as well, else print 1
               if cheaters are the first to die, print 3, and then if everyone dies, print 0 as well, else print 1
               if everyone coexists, then print 5
            */

            let outcome;

            if (total_broadcaster_system === 0 || total_tetherer_system === 0 || total_cheater_system === 0) {
                if (!extinction_event) {
                    extinction_event = true;
                    extinction_time = sim.time
                }
                if (total_broadcaster_system === 0) {
                    outcome = 1;
                }
                else if (total_tetherer_system === 0) {
                    outcome = 2;
                }
                else if (total_cheater_system === 0) {
                    outcome = 3;
                }
            }

            if (extinction_event && sim.time > extinction_time + 20000) {
                let final_outcome = 1;
                if (total_broadcaster_system + total_tetherer_system + total_cheater_system === 0) {
                    final_outcome = 0
                }   
                let out_1 = destruction_counters + '\n'
                fs.appendFileSync(destruction_count, out_1)
                const resultObject = {
                    run_num: run_num,
                    outcome: outcome,
                    final_outcome: final_outcome,
                    starting_polymer_value: starting_polymer_value,
                    migration_rate: migration,
                    destruction_rate: destruction
                };
                
                fs.appendFileSync(final_result,'[' + JSON.stringify(resultObject) + ']' + '\n');
                
                process.exit();
            } 

            else if (sim.time === config.maxtime - 1) {
                const resultObject = {
                    run_num: run_num,
                    outcome: 5,
                    starting_polymer_value: starting_polymer_value,
                    migration_rate: migration,
                    destruction_rate: destruction
                };
            
                fs.appendFileSync(final_result,'[' + JSON.stringify(resultObject) + ']' + '\n');
            } 
        }

        // now calling cacatoo's update() functions will lead to them updating in the shuffled order

        colonies.update = function() {
            //console.log(x);  
            colonies.asynchronous();       // asynchronously update the gridpoints in a random order. this introduces stochasticity and ensures the gridpoints are not updated in a 'set' order
            if (sim.time % (diff_interval_final) === 0 && (mar_diffusion === true)) {
                (colonies.MargolusDiffusion())   // implement Margolus diffusion
            }
            else if (sim.time % (diff_interval_final) === 0 && (mix === true)) {
                (colonies.perfectMix())     // this is pretty much useless here since we do not want individual grids to have well-mixed colonies
            }

            //colonies.plotPopsizes('type',['b', 'c'])
            //console.log("bleh", x)  

            let sum_p = 0;
            let sum_h = 0;
            let sum_c = 0;
            let sum_internal_resources_p = 0;
            let sum_internal_resources_h = 0;
            let sum_internal_resources_c = 0;
            let avg_int_p = 0;
            let avg_int_h = 0;
            let avg_int_c = 0;
            let monomer_per_gp = 0.0;
            let polymer_per_gp = 0.0;

            for (let i = 0; i < colonies.nc; i++) {
                for (let j = 0; j < colonies.nr; j++) {
                    if (colonies.grid[i][j].type === 'b') {
                        sum_p++;
                        sum_internal_resources_p += colonies.grid[i][j].internal_resources;
                    }
                    else if (colonies.grid[i][j].type === 't') {
                        sum_h++;
                        sum_internal_resources_h += colonies.grid[i][j].internal_resources;
                    }
                    else if (colonies.grid[i][j].type === 'c') {
                        sum_c++;
                        sum_internal_resources_c += colonies.grid[i][j].internal_resources;
                    }
                }
            }

            for (let x = 0; x < resources.nc; x++) {
                for (let y = 0; y < resources.nr; y++) {
                    monomer_per_gp += resources.grid[x][y].monomer_count      // tracks the monomer units present per grid point
                    polymer_per_gp += resources.grid[x][y].polymer_count
                }
            }  

            avg_int_p = (sum_p === 0) ? 0 : sum_internal_resources_p / sum_p;
            avg_int_h = (sum_h === 0) ? 0 : sum_internal_resources_h / sum_h;
            avg_int_c = (sum_c === 0) ? 0 : sum_internal_resources_c / sum_c;

            broadcaster_populations[x] = sum_p;
            tetherer_populations[x] = sum_h;
            cheater_populations[x] = sum_c;
            broadcaster_internal_resources[x] = sum_internal_resources_p;
            tetherer_internal_resources[x] = sum_internal_resources_h;
            cheater_internal_resources[x] = sum_internal_resources_c;

            // the file writing format is - timestamp / grid uID / braodcaster pop size / broadcaster avg biomass / tetherer pop size / tetherer avg biomass / avg polymer / avg monomer

            if (sim.time % data_interval === 0) {
                let out_1 = sim.time + '\t\t'
                let out_2 = resources.uID + '\t\t'
                let out_3 = sum_p + '\t\t' 
                let out_4 = avg_int_p.toExponential(2) + '\t\t'
                let out_5 = sum_h + '\t\t'
                let out_6 = avg_int_h.toExponential(2) + '\t\t'
                let out_7 = sum_c + '\t\t'
                let out_8 = avg_int_c.toExponential(2) + '\t\t'
                let out_9 = (polymer_per_gp/(resources.nc*resources.nr)).toExponential(2) + '\t\t'
                let out_10 = (monomer_per_gp/(resources.nc*resources.nr)).toExponential(2) + '\n'
                fs.appendFileSync(pop_output_files[x], out_1)
                fs.appendFileSync(pop_output_files[x], out_2)
                fs.appendFileSync(pop_output_files[x], out_3)
                fs.appendFileSync(pop_output_files[x], out_4)
                fs.appendFileSync(pop_output_files[x], out_5)
                fs.appendFileSync(pop_output_files[x], out_6)
                fs.appendFileSync(pop_output_files[x], out_7)
                fs.appendFileSync(pop_output_files[x], out_8)
                fs.appendFileSync(pop_output_files[x], out_9)
                fs.appendFileSync(pop_output_files[x], out_10)
            }

            if (run_num === 1) {
                for (let timestamp in timestamps_of_interest) {
                    if (sim.time === timestamps_of_interest[timestamp] && grid_capture === false) {
                        grid_capture = true;
                        time_at_grid_capture = sim.time
                        //console.log(time_at_grid_capture)
                        grid_index_for_images = x
                    }
                }

                if (grid_capture === true) {
                    if (sim.time === (time_at_grid_capture + image_interval)) {
                        grid_capture = false;
                        //console.log("done")
                        grid_index_for_images = undefined;
                    };    

                    //console.log(x, grid_index_for_images)    
                    if (take_png === true && x === grid_index_for_images && grid_capture === true) {
                        //console.log(grid_index_for_images)
                        //console.log(x, grid_index_for_images+'lol')                    
                        for (let i = 0; i < colonyGridModels[grid_index_for_images].nc; i++) {
                            for (let j = 0; j < colonyGridModels[grid_index_for_images].nr; j++) {
                                let cell = colonyGridModels[grid_index_for_images].grid[i][j]
                                //console.log(i, j, cell.type);
                                //console.log(`Processing cell (${i},${j}): type=${cell.type}`);
                                color = '#000000'
                                if (cell.type === 'b') {
                                    color = '#9F2B68';  // amaranth
                                } 
                                else if (cell.type === 't') {
                                    color = '#FFA700';  // chrome yellow
                                }
                                else if (cell.type === 'c') {
                                    color = '#00FFFF'; // cyan
                                }
                                //console.log(`Using fill color: ${color}`)
                        
                                ctx.fillStyle = color
                                ctx.fillRect(i*cellSize, j*cellSize, cellSize, cellSize)
                            }
                        }
                        //console.log(grid_capture)
                        b = canvas.toBuffer('image/png')
                        //console.log(grid_index_for_images)
                        writeFileSync(join(movie, `grid_${resourcesGridModels[grid_index_for_images].uID}_${sim.time}.png`), b)
                        ctx.clearRect(0, 0, canvasWidth, canvasHeight);    // clear the canvas after each image genration to free up RAM
                        //console.log(grid_index_for_images+"lmao")
                        b = null;  // to mark the buffer for garbage collection
                    }
                    else if (take_txt === true && x === grid_index_for_images && grid_capture === true) {
                        let heatmap_file = heatmaps+"/heatmap_type"+sim.time+"_Grid"+grid_index_for_images+".txt"
                        for (let i = 0; i < colonyGridModels[grid_index_for_images].nc; i++) {
                            for (let j = 0; j < colonyGridModels[grid_index_for_images].nr; j++) {
                                let cell = colonyGridModels[grid_index_for_images].grid[i][j]
                                let point = i+'\t'+j+'\t'+cell.type+'\n'
                                fs.appendFileSync(heatmap_file,point)
                            }
                        }
                    }
                }    
            }
        }
        
        resources.update = function() {
            resources.asynchronous();
            //console.log(x, "meh");

            // destorying a grid. this is implementing the Intermediate Disturbance Hypothesis (in a way, as what is a 'disturbance' is subjective). each individual project's destruction_rate should be tested
            // to find the right frequency of disturbances and check if coexistence occurs or not

            if (this.rng.genrand_real1() < destruction_rate) {
                let destroyed_grid = {                    // create an object that stores the timestamp of destruction and the grid model
                    colony: colonies,
                    resource: resources,
                    destroy_time: sim.time
                }
                // iterate over every grid point and destroy!
                for (let x = 0; x < colonies.nc; x++) {
                    for (let y = 0; y < colonies.nr; y++) {
                        // continue your mayhem
                        colonies.grid[x][y].type = 0.0;
                        colonies.grid[x][y].Km_polymer = 0.0;
                        colonies.grid[x][y].Vmax_polymer = 0.0;
                        colonies.grid[x][y].Km_monomer = 0.0;
                        colonies.grid[x][y].Vmax_monomer = 0.0;
                        colonies.grid[x][y].internal_resources = 0.0;
                    }
                }

                /*let monomer_per_gp_at_destruction = 0.0
                //let polymer_per_gp_at_destruction = 0.0 

                for (let x = 0; x < resources.nc; x++) {
                    for (let y = 0; y < resources.nr; y++) {
                        monomer_per_gp_at_destruction += resources.grid[x][y].monomer_count      // tracks the monomer units present per grid point
                        //polymer_per_gp_at_destruction += resources.grid[x][y].polymer_count
                    }
                }*/    
                let out_1 = sim.time + '\t\t';
                let out_2 = `grid ${x}` + '\n';
                fs.appendFileSync(destroyed_grids_file, out_1);
                fs.appendFileSync(destroyed_grids_file, out_2);
                destroyedGrids.push(destroyed_grid)  // store the destroyed grid in an array to be used later for refreshing
                // finally clear the resources as well. don't worry, the refresh logic ensures a fresh, uncolonised particle is created after destruction after a user-defined time gap
                
                let out_4 = 'destroyed' + '\n'
                //let out_5 = (monomer_per_gp_at_destruction/(resources.nc*resources.nr)).toExponential(2) + '\n'
                fs.appendFileSync(pop_output_files[x], out_1)
                fs.appendFileSync(pop_output_files[x], out_4)
                //fs.appendFileSync(pop_output_files[x], out_5)
                destruction_counter++;
                let out_3 = 'destroyed' + '\n'
                fs.appendFileSync(migrants_output_files[x], out_3)
                migration_counter = 0;
                invasion_counter = 0;
                destruction_counters[x] = destruction_counter;

                sim.initialGrid(resources, 'polymer_count', 0.0, 1.0);
                sim.initialGrid(resources, 'monomer_count', 0.0, 1.0);
            }
            
            refreshDestroyedGrids();

            // one of the most important code blocks here. migration is an important aspect of movement ecology especially when it comes to metapopulation dynamics
                      
            if (sim.time % migration_time === 0) {        // choose to perform migration after a set number of time steps or every time step
                let destination_index;               
                let currentGrid = colonies;
                let p_migrants = 0;
                let h_migrants = 0;
                let c_migrants = 0;
                let p_migrants_internal_resources = 0;
                let h_migrants_internal_resources = 0;
                let c_migrants_internal_resources = 0;

                for (let i = 0; i < currentGrid.nc; i++) {
                    for (let j = 0; j < currentGrid.nr; j++) {
                        if (currentGrid.grid[i][j].type !== 0) {
                            if (this.rng.genrand_real1() < migration_rate) {   // check if a grid point has a bacterium and then check if it will migrate
                                do {
                                    destination_index = Math.floor(Math.random() * colonyGridModels.length);
                                } while (destination_index === x);             // keep generating a new value while the destination index is the 'home' grid
                                let destinationGrid = colonyGridModels[destination_index];
                                migration_counter++;

                                // copy all the properties of the migrant
                                let migrant_cell = currentGrid.grid[i][j].type;
                                let migrant_vmax_polymer = currentGrid.grid[i][j].Vmax_polymer
                                let migrant_km_polymer = currentGrid.grid[i][j].Km_polymer
                                let migrant_vmax_monomer = currentGrid.grid[i][j].Vmax_monomer
                                let migrant_km_monomer = currentGrid.grid[i][j].Km_monomer
                                let migrant_fitness = currentGrid.grid[i][j].internal_resources

                                // find a destination grid point. does not matter if it is occupied. we are allowing invasion. this can be changed as per the user's requirements
                                let dest_i = Math.floor(Math.random() * destinationGrid.nc);
                                let dest_j = Math.floor(Math.random() * destinationGrid.nr);

                                // write the migrant's information on to a file. this writes its fitness, the destination and source grids and whether it invaded or not
                                // alternatively, since the system as a whole is well-mixed, tracking the source and destination grids becomes useless as there is no distinct to and fro 
                                // pattern emerging here.

                                if (currentGrid.grid[i][j].type === 'b') {
                                    p_migrants++;
                                    p_migrants_internal_resources += currentGrid.grid[i][j].internal_resources;
                                }
                                else if (currentGrid.grid[i][j].type === 't') {
                                    h_migrants++;
                                    h_migrants_internal_resources += currentGrid.grid[i][j].internal_resources;
                                }
                                else if (currentGrid.grid[i][j].type === 'c') {
                                    c_migrants++;
                                    c_migrants_internal_resources += currentGrid.grid[i][j].internal_resources;
                                }
                                
                                if (destinationGrid.grid[dest_i][dest_j].type !== 0) {
                                    invasion_counter++;
                                }
           
                                // time to move the migrant cell to the destination cell
                                destinationGrid.grid[dest_i][dest_j].type = migrant_cell;
                                destinationGrid.grid[dest_i][dest_j].Vmax_polymer = migrant_vmax_polymer;
                                destinationGrid.grid[dest_i][dest_j].Km_polymer = migrant_km_polymer;
                                destinationGrid.grid[dest_i][dest_j].Vmax_monomer = migrant_vmax_monomer;
                                destinationGrid.grid[dest_i][dest_j].Km_monomer = migrant_km_monomer;
                                destinationGrid.grid[dest_i][dest_j].internal_resources = migrant_fitness  - migration_cost;

                                // and clear the grid point at the source grid
                                currentGrid.grid[i][j].type = 0;
                                currentGrid.grid[i][j].Vmax_polymer = 0;
                                currentGrid.grid[i][j].Km_polymer = 0;
                                currentGrid.grid[i][j].Vmax_monomer = 0;
                                currentGrid.grid[i][j].Km_monomer = 0;
                                currentGrid.grid[i][j].internal_resources = 0;
                            }                            
                        }
                    }
                }
                
                let out_1 = sim.time + '\t'
                let out_2 = p_migrants + '\t'
                let out_3 = h_migrants + '\t'
                let out_4 = c_migrants + '\t'
                let out_5 = (p_migrants_internal_resources/p_migrants).toExponential(2) + '\t'
                if (isNaN(out_4)) {out_5 = 0 + '\t'}
                let out_6 = (h_migrants_internal_resources/h_migrants).toExponential(2) + '\t'
                if (isNaN(out_5)) {out_6 = 0 + '\t'}
                let out_7 = (c_migrants_internal_resources/c_migrants).toExponential(2) + '\t'
                if (isNaN(out_7)) {out_7 = 0 + '\t'}
                let out_8 = invasion_counter + '\t'
                let out_9 = migration_counter + '\n'
                fs.appendFileSync(migrants_output_files[x], out_1)
                fs.appendFileSync(migrants_output_files[x], out_2)
                fs.appendFileSync(migrants_output_files[x], out_4)
                fs.appendFileSync(migrants_output_files[x], out_3)
                fs.appendFileSync(migrants_output_files[x], out_5)
                fs.appendFileSync(migrants_output_files[x], out_6)
                fs.appendFileSync(migrants_output_files[x], out_7)
                fs.appendFileSync(migrants_output_files[x], out_8)
                fs.appendFileSync(migrants_output_files[x], out_9)
                migration_counter = 0;
                invasion_counter = 0;               
            }

            if (sim.time === config.maxtime - 1) {
                let out_1 = x + '\t'
                let out_2 = destruction_counter + '\n'
                fs.appendFileSync(destruction_count, out_1)
                fs.appendFileSync(destruction_count, out_2)
            }
        }
    }
}
//console.log(grid_index)

// call on your functions

defineNextStates();
defineUpdates();

// time to tango

//sim.pause = true;
//sim.addButton("Play/Pause", function () { sim.toggle_play() })
sim.start()
