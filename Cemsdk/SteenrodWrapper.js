'use strict';
importScripts("CSteenrod.js");
// let ctypes = require("ctypes");

function addClassCallback(hom_deg, int_deg, cocycle) { 
    self.postMessage({"cmd" : "addClass", "x" : int_deg - hom_deg, "y": hom_deg});
    // display.updateBatch();
}

function addStructlineCallback(
    source_hom_deg, source_int_deg, source_idx, 
    target_hom_deg, target_int_deg, target_idx
){
    self.postMessage({"cmd" : "addStructline", 
        "source" : {"x" : source_int_deg - source_hom_deg, "y": source_hom_deg, "idx": source_idx},
        "target" : {"x" : target_int_deg - target_hom_deg, "y": target_hom_deg, "idx": target_idx}
    })
}


let addClassCallbackPtr = Module.addFunction(addClassCallback);
let addStructlineCallbackPtr = Module.addFunction(addStructlineCallback);
let cdoRes = cwrap("doResolution", 'void', ['number', 'pointer', 'pointer'])

function doRes(degree){
    cdoRes(degree, addClassCallbackPtr, addStructlineCallbackPtr);
}

let runtimePromise = new Promise(function(resolve, reject){
    Module.onRuntimeInitialized = function(){
        resolve();
    }; 
});

self.onmessage = function(msg){
    runtimePromise.then(() => doRes(msg.data.degree));
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