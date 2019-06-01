'use strict';
importScripts("CSteenrod.js");
importScripts("CSteenrodWrappers.js");

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





function constructModule(module, cAlgebra){
    console.log(module);
    let gensArray = new Uint32Array(Object.values(module.gens));
    let max_degree = Math.max(...Object.values(module.gens));
    let graded_dimension = new Uint32Array(max_degree + 1);
    let basis_element_indices = {};
    let index_to_basis_element = {};
    for(let [b, degree] of Object.entries(module.gens)){
        let index = graded_dimension[degree];
        graded_dimension[degree] ++;
        basis_element_indices[b] = index;
        index_to_basis_element[degree] = index_to_basis_element[degree] || {}
        index_to_basis_element[degree][index] = b;
    }
    let c_array_offset = Module._malloc(4 * Math.max((max_degree + 1), ...graded_dimension));
    Module.HEAPU32.set(graded_dimension, c_array_offset/4);
    let cModule = cFiniteDimensionalModule_construct(cAlgebra, max_degree, c_array_offset);
        for(let {op, input, output} of module.milnor_actions){
        op = op.concat([0,0,0]);
        let op_degree = op[0] + 3*op[1] + 7*op[2]; // TODO: fix me.
        let input_degree = module.gens[input];
        let input_index = basis_element_indices[input];
        let output_degree = input_degree + op_degree;        
        let op_index = 0; //cMilnorBasisElement.toIndex(algebra, op);        
        let output_vector = new Uint32Array(graded_dimension[output_degree]);
        console.log(output);
        for( let {gen, coeff} of output) {
            console.log(gen, coeff);
            output_vector[basis_element_indices[gen]] = coeff
        }
        console.log(output_vector);
        Module.HEAPU32.set(output_vector, c_array_offset/4);
        cFiniteDimensionalModule_setAction(
            cModule, 
            op_degree, op_index,
            input_degree, input_index,
            c_array_offset
        )
    }
    Module._free(c_array_offset);
    console.log("Finished module!");
    return cModule;
}


function doRes(degree){
    cdoRes(degree, addClassCallbackPtr, addStructlineCallbackPtr);
}

let runtimePromise = new Promise(function(resolve, reject){
    Module.onRuntimeInitialized = function(){
        resolve();
    }; 
});

self.onmessage = function(msg){
    console.log(msg);
    runtimePromise.then(() => {
        let p = msg.data.module.p;
        console.log(msg.data.module);
        let generic = msg.data.module.generic;
        let max_degree = msg.data.max_degree;
        if(typeof(p) != 'number'){
            console.log("p not a number, quitting.");
            return;
        }        
        if(typeof(max_degree) != 'number'){
            console.log("max_degree not a number, quitting.");
            return;
        }
        cinitializePrime(p);
        let cAlgebra = cMilnorAlgebra_construct(p, generic, null);    
        cMilnorAlgebra_generateBasis(cAlgebra, max_degree);
        let cModule = constructModule(msg.data.module, cAlgebra);
        let cResolution = cResolution_construct(cModule, max_degree, addClassCallbackPtr, addStructlineCallbackPtr);
        cResolution_resolveThroughDegree(cResolution, max_degree);
    });
};





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