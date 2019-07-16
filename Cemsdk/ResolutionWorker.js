'use strict';
importScripts("CSteenrod.js");
importScripts("CSteenrodWrappers.js");

function constructCString(string){
    let offset = Module._malloc(string.length + 1);
    javascriptStringToC(offset, string);
    return offset;
}


function _arrayBufferToBase64( buffer ) {
    var binary = '';
    var bytes = new Uint8Array( buffer );
    var len = bytes.byteLength;
    for (var i = 0; i < len; i++) {
        binary += String.fromCharCode( bytes[ i ] );
    }
    return self.btoa( binary );
}

function javascriptStringToC(offset, string){
    for(let i = 0; i < string.length; i++){
        // Hopefully all of the characters are <= 255
        Module.HEAPU8[offset + i] = string.charCodeAt(i);
    }
    Module.HEAPU8[offset + string.length + 1] = 0;
}

function cStringToJavascript(offset, length){
    let result = [];
    for(let i = 0; i < length; i++){
        result.push(Module.HEAPU8[offset + i]);
    }
    return String.fromCharCode(...result);
}

function FreeModuleElement_to_display_string(free_module, degree, vector) {
    let cResult_json_offset = Module._malloc(2000);
    let len = cFreeModule_element_toJSONString(cResult_json_offset, free_module, degree, vector);
    let s = "";
    for (let i = 0; i < len; ++i){
        s += String.fromCharCode(Module.HEAPU8[cResult_json_offset+i]);
    } 
    let result = FreeModuleElement_json_to_display_string(s);
    Module._free(cResult_json_offset);
    return result;
}

function FreeModuleElement_json_to_display_string(json){
    let result = JSON.parse(json).map((entry) => {
        let output = "";
        if(entry.coeff != 1){
            output = `${entry.coeff}*`;
        }
        let op_string = entry.op_deg > 0 ? `${entry.op_str}*`  : "";
        output = `${output}${op_string}x_{${entry.gen_deg},${entry.gen_idx}}`
        return output;
    })
    if(result.length == 0){
        result = "0";
    } else {
        result = result.join(" + ");
    }
    return result;
}

const sizeof_uint = Uint32Array.BYTES_PER_ELEMENT;

let t0 = performance.now();
let t_last = t0;
function getTime(){
    let t_cur = performance.now();
    let duration = (t_cur - t_last) / 1000;
    t_last = t_cur;
    return duration;
}
function getTotalTime(){
    let t_cur = performance.now();
    return (t_cur - t0) / 1000;
}
let last_timestamp_stem = 0;
let timestamp_interval = 10;
function getCallbacks(){
    function addClassCallback(hom_deg, int_deg, cocycle) { 
        if(int_deg >= last_timestamp_stem + timestamp_interval){
            let total_time = getTotalTime().toFixed(2) + " seconds";
            let last_ten_stems_time = getTime().toFixed(2) + " seconds";
            console.log(`Time to compute stems ${last_timestamp_stem} to ${last_timestamp_stem + timestamp_interval}: ${last_ten_stems_time}`);
            console.log(`Total time to compute first ${int_deg} stems: ${total_time}`);
            last_timestamp_stem += timestamp_interval;
        }
        self.postMessage({"cmd" : "addClass", "x" : int_deg - hom_deg, "y": hom_deg});
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
    
    // TODO: remove these addFunctions as part of cleanup?
    let addClassCallbackPtr = Module.addFunction(addClassCallback);
    let addStructlineCallbackPtr = Module.addFunction(addStructlineCallback);
    return {addClassPtr : addClassCallbackPtr, addStructlinePtr : addStructlineCallbackPtr};
}



function constructAlgebra(algebraData){
    let p = algebraData.p;
    let generic = algebraData.generic;
    cinitializePrime(p);
    // Module.print(algebraData);
    let cAlgebra;
    if(!algebraData.algebra || algebraData.algebra.type.toLowerCase() == "milnor"){
        let cProfile = 0;
        if(algebraData.algebra && algebraData.algebra.profile){
            let profile = algebraData.algebra.profile;
            let q_part = profile.q_part || [];
            let p_part = profile.p_part;
            let truncated = profile.truncated;        
            let c_qpart_offset = Module._malloc(sizeof_uint * (q_part.length + p_part.length));
            let c_ppart_offset = c_qpart_offset + sizeof_uint*q_part.length;
            Module.HEAPU32.set(new Uint32Array(q_part), c_qpart_offset/sizeof_uint);
            Module.HEAPU32.set(new Uint32Array(p_part), c_ppart_offset/sizeof_uint);
            cProfile = cProfile_construct(p != 2, q_part.length, c_qpart_offset, p_part.length, c_ppart_offset, truncated);
            Module._free(c_qpart_offset);
        }
        cAlgebra = cMilnorAlgebra_construct(p, generic, cProfile);
        cProfile_free(cProfile); 
        cMilnorAlgebra_generateBasis(cAlgebra, algebraData.max_degree);
    } else {
        let unstable = algebraData.algebra.unstable || false;
        cAlgebra = cAdemAlgebra_construct(p, generic, unstable);
        cAdemAlgebra_generateBasis(cAlgebra, algebraData.max_degree);    
    }
    return cAlgebra;
}


function constructFiniteDimensionalModule(module, cAlgebra){
    let p = module.p;
    let min_basis_degree = Math.min(...Object.values(module.gens))
    let max_basis_degree = Math.max(...Object.values(module.gens)) - min_basis_degree + 1;
    let graded_dimension = new Uint32Array(max_basis_degree);
    let basis_element_indices = {};
    let index_to_basis_element = {};
    for(let [b, degree] of Object.entries(module.gens)){
        let index = graded_dimension[degree - min_basis_degree];
        graded_dimension[degree - min_basis_degree] ++;
        basis_element_indices[b] = index;
        index_to_basis_element[degree - min_basis_degree] = index_to_basis_element[degree - min_basis_degree] || {}
        index_to_basis_element[degree - min_basis_degree][index] = b;
    }
    let c_module_name = constructCString(module.file_name);
    let c_array_offset = Module._malloc(sizeof_uint * Math.max(max_basis_degree, ...graded_dimension));
    Module.HEAPU32.set(graded_dimension, c_array_offset/sizeof_uint);
    let cModule = cFiniteDimensionalModule_construct(cAlgebra, c_module_name, min_basis_degree, min_basis_degree + max_basis_degree, c_array_offset);
    Module._free(c_module_name);
    for(let {op, input, output} of module.milnor_actions){
        let op_degree;
        if(p==2){
            op = op.concat([0,0,0]);
            op_degree = op[0] + 3*op[1] + 7*op[2]; // TODO: fix me.
        } else {
            let opQ = op[0];
            let opP = op[1].concat([0,0]);
            op_degree = 2*(p-1)*opP[0] + opQ.length;
        }
        if(op_degree == 0){
            continue;
        }
        let input_degree = module.gens[input];
        let input_index = basis_element_indices[input];
        let output_degree = input_degree + op_degree;        
        let op_index = 0; //cMilnorBasisElement.toIndex(algebra, op);        
        let output_vector = new Uint32Array(graded_dimension[output_degree - min_basis_degree]);
        for( let {gen, coeff} of output) {
            output_vector[basis_element_indices[gen]] = coeff
        }
        Module.HEAPU32.set(output_vector, c_array_offset/sizeof_uint);

        cFiniteDimensionalModule_setAction(
            cModule, 
            op_degree, op_index,
            input_degree, input_index,
            c_array_offset
        )
    }
    Module._free(c_array_offset);
    let result = {cModule : cModule};
    return result;
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
    runtimePromise.then(() => {
        message_handlers[msg.data.cmd](msg.data);
    });
};

let message_handlers = {};

message_handlers["resolve"] = function resolve(data){
    let p = data.module.p;
    let max_degree = data.max_degree;
    max_degree++;
    if(typeof(p) != 'number'){
        console.log("p not a number, quitting.");
        return;
    }        
    if(typeof(max_degree) != 'number'){
        console.log("max_degree not a number, quitting.");
        return;
    }
    
    let algebraData = data.algebra || {};
    let moduleData = data.module;
    algebraData.p = moduleData.p;
    algebraData.generic = moduleData.generic;
    algebraData.max_degree = max_degree - Math.min(...Object.values(moduleData.gens));
    moduleData.max_degree = max_degree;
    algebraData.algebra = moduleData.algebra;
    let cAlgebra = constructAlgebra(algebraData);
    let module = constructFiniteDimensionalModule(data.module, cAlgebra);
    let callbacks = getCallbacks();
    t0 = performance.now();
    self.p = p;
    self.cResolution = cResolution_construct(module.cModule, max_degree, callbacks.addClassPtr, callbacks.addStructlinePtr);
    self.cResWithMaps = cResolutionWithMapsToUnitResolution_construct(self.cResolution, self.cResolution, 2);
    cresolveThroughDegree(self.cResWithMaps, max_degree);

    // self.x_homological_degree = 1;
    // self.x_internal_degree = 2;
    // self.x_idx = 0;
    // self.f = cResolutionHomomorphism_construct(self.cResolution, self.cResolution, x_homological_degree, x_internal_degree);
    // let v = cVector_construct(p, 1, 0);
    // cVector_setEntry(v, 0, 1);
    // cResolutionHomomorphism_setBaseMap(self.f, x_internal_degree, x_idx, v);
    // cResolutionHomomorphism_baseMapReady(self.f, 1000);
    // cResolutionHomomorphism_extend(f, max_degree - 2, max_degree - 2);
    // cVector_free(v); v = null;
    // message_handlers["serialize"](0);
};

message_handlers["get_cocycle"] = function getCocycle(data){
    let homological_degree = data.y;
    let degree = data.x + data.y;
    let index = data.idx;
    let cf = cResolution_getDifferential(self.cResolution, homological_degree);
    let cTarget = cFreeModuleHomomorphism_getTarget(cf);
    let dimension = cModule_getDimension(cTarget, degree);
    let cResult_vector = cVector_construct(self.p, dimension, 0);    
    cFreeModuleHomomorphism_applyToGenerator(cf, cResult_vector,  1, degree, index);
    cVector_print("%s\n", cResult_vector);
    let result = FreeModuleElement_to_display_string(cTarget, degree, cResult_vector)
    self.postMessage({
        "cmd" : "cocycle_result",
        "value" : result
    });
};

message_handlers["get_eta_map"] = function get_eta_map(data){
    let homological_degree = data.y;
    let degree = data.x + data.y;
    let output_degree = degree - self.x_internal_degree;
    let index = data.idx;
    let cf = cResolutionHomomorphism_getMap(self.f, homological_degree);
    let cTarget = cFreeModuleHomomorphism_getTarget(cf);
    let dimension = cModule_getDimension(cTarget, output_degree);
    let cResult_vector = cVector_construct(self.p, dimension, 0);    
    cFreeModuleHomomorphism_applyToGenerator(cf, cResult_vector,  1, degree, index);
    cVector_print("%s\n", cResult_vector);
    let result = FreeModuleElement_to_display_string(cTarget, output_degree, cResult_vector)
    self.postMessage({
        "cmd" : "cocycle_result",
        "value" : result
    });
};


message_handlers["serialize"] = function serialize(data){
    let serialized_resolution = cResolution_serialize(self.cResolution);
    let serialized_resolution_json_data = cSerializedResolution_getJSONData(serialized_resolution);
    let serialized_resolution_json_length = cSerializedResolution_getJSONSize(serialized_resolution);
    let serialized_resolution_binary_data = cSerializedResolution_getBinaryData(serialized_resolution);
    let serialized_resolution_binary_length = cSerializedResolution_getBinarySize(serialized_resolution);
    let json_string = cStringToJavascript(serialized_resolution_json_data, serialized_resolution_json_length);
    let buffer = new ArrayBuffer(serialized_resolution_binary_length);
    for(let i=0; i<serialized_resolution_binary_length; i++){
        buffer[i] = Module.HEAPU8[serialized_resolution_binary_data + i];
    }
    // Module.print("byte_length: ", buffer.byteLength);
    self.postMessage({
        "cmd" : "serialized_matrices",
        "binary_data" : buffer,
        "json_data" : json_string
    }, [buffer]);
    // Module.print("byte_length: ", buffer.byteLength);
    // Module.print(json_string);
}


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