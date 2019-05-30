'use strict';
// let ctypes = require("ctypes");
let sseq = new Sseq();
sseq.display();

function addClassCallback(hom_deg, int_deg, cocycle) { 
    let c = sseq.addClass(int_deg - hom_deg, hom_deg);
    classesTable[int_deg - hom_deg][hom_deg].push(c);
    sseq.update();
}

function addStructlineCallback(
    source_hom_deg, source_int_deg, source_idx, 
    target_hom_deg, target_int_deg, target_idx
){
    let source = classesTable[source_int_deg - source_hom_deg][source_hom_deg][source_idx];
    let target = classesTable[target_int_deg - target_hom_deg][target_hom_deg][target_idx];
    sseq.addStructline(source, target);
}


let addClassCallbackPtr = Module.addFunction(addClassCallback);
let addStructlineCallbackPtr = Module.addFunction(addStructlineCallback);
let cdoRes = cwrap("doResolution", 'void', ['number', 'pointer', 'pointer'])

function doRes(degree){
    window.classesTable = Array(degree).fill(0).map(x => Array(degree).fill(0).map(x => []));
    cdoRes(degree, addClassCallbackPtr, addStructlineCallbackPtr);
}


Module.onRuntimeInitialized = function(){
    doRes(50);
}; 
// sleep(1000).then( () => doRes(20));




// // typedef struct Resolution_s {
// //     Algebra *algebra; // 0
// //     Module *module;   // 4
// //     void (*callback)(struct Resolution_s* res, uint hom_deg, uint int_deg); // 8
// //     uint max_degree; // 12
// //     FreeModule **modules; // The index into resolution_modules is homological_degree + 1. // 16
// //     FreeModuleHomomorphism **differentials;// 20
// //     int *internal_degree_to_resolution_stage;// 24
// // } Resolution;
// let resolutionStruct = {};
// resolutionStruct.algebraOffset = 0;
// resolutionStruct.moduleOffset = 8;
// resolutionStruct.moduleOffset = 16;
// resolutionStruct.moduleOffset = 20;

// Module.print( `Res ${res} finished (${hom_deg}, ${int_deg}).`
//     + `There ${num_gens != 1 ? "are" : "is"} ${num_gens} homology generator${num_gens != 1 ? "s" : ""}.`);
// // for(let offset = 0; offset < 80; offset += 4){
// //     Module.print(`res + ${offset}:`, Module.getValue(res+offset));
// // }
// // let max_degree = Module.getValue(res + 12);
// // Module.print(`max_degree: ${max_degree}`);
// // // maybe this is internal_degree_to_resolution_stage?
// // let internal_degree_to_resolution_stage = Module.getValue(res+24);
// // for(let index = 0; index < 5 ; index ++){
// //     Module.print(
// //         `int_deg_to_res_stage${index}: `  
// //         + Module.getValue(internal_degree_to_resolution_stage + index*4).toString()
// //     );
// // }
// // Module.print("\n\n\n");

// for(let i = 0; i < num_gens; i++){