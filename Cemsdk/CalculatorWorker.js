'use strict';
importScripts("CSteenrod.js");
importScripts("CSteenrodWrappers.js");
const sizeof_uint = 4;
const string_buffer_size = 10000;
let algebra; 
let degree_caps = {
    2 : 105,
    3 : 313,
    5 : 882,
    7 : 1537
}


function constructAlgebra(algebraData){
    // if(algebra){
    //     Module.print("Already constructed algebra...");
    //     return;
    // }
    let p = algebraData.p;
    let generic = algebraData.generic;
    let q = generic ? 2*p-2 : 1;
    let max_degree = 20*q + 1;
    let unstable = true;
    cinitializePrime(p);
    let c_algebra;
    if(algebraData.type.toLowerCase() === "milnor"){
        let cProfile = 0;
        c_algebra = cMilnorAlgebra_construct(p, generic, cProfile);
    } else if(algebraData.type.toLowerCase() === "adem"){
        c_algebra = cAdemAlgebra_construct(p, generic, unstable);
    } else {
        Module.print("Unrecognized algebra type!");
        return undefined;
    }
    let algebra = {};
    algebra.p = p;
    algebra.generic = generic;
    algebra.max_degree = 0;
    algebra.type = algebraData.type
    algebra.c_algebra = c_algebra;
    algebra.string_buffer = Module._malloc(string_buffer_size);
    generateBasis(algebra, max_degree);
    return algebra;
}

function generateBasis(algebra, max_degree){
    if(max_degree < algebra.max_degree){
        return;
    }
    self.postMessage({"cmd" : "adem_extending_basis", "old_max" : algebra.max_degree, "new_max" : max_degree - 1});
    algebra.max_degree = max_degree;
    if(algebra.type.toLowerCase() === "milnor"){
        cMilnorAlgebra_generateBasis(algebra.c_algebra, max_degree);    
    } else if(algebra.type.toLowerCase() === "adem"){
        cAdemAlgebra_generateBasis(algebra.c_algebra, max_degree);    
    } else {
        Module.print("Unrecognized algebra type!");
        return undefined;
    }
    self.postMessage({"cmd" : "adem_extending_basis_done"});
}

let runtimePromise = new Promise(function(resolve, reject){
    Module.onRuntimeInitialized = function(){
        resolve();
    }; 
});

self.onmessage = function(msg){
    runtimePromise.then(() => {
        handlers[msg.data.cmd](msg.data);
    });
};

let handlers = {};

function ademSequenceTojsAdemMonomial(q, input){
    input = input.trim();
    let Ps = [];
    let bocksteins = 0;
    let degree = 0;
    for(let s of input.split(" ").values()){
        if(s === "b"){
            degree ++;
            bocksteins |= (1 << Ps.length);
        } else if(s === ""){
            continue;
        } else if(isNaN(s)){
            return { "error" : s, "position" : input.indexOf(s)};
        } else {
            s = parseInt(s);
            degree += s * q;
            Ps.push(s);
        }
    }
    return {
        "degree" : degree,
        "bocksteins" : bocksteins,
        "Ps" : Ps
    };
}

function jsAdemMonomialToOutput(Ps, bocksteins, p){
    let simple_result = [];
    let latex_result = []
    let P_or_Sq = p === 2 ? "Sq" : "P";
    for(let [idx,P] of Ps.entries()){
        if((bocksteins >> idx) & 1){
            simple_result.push("b");
            latex_result.push("\\beta");
        }
        simple_result.push(P.toString());
        latex_result.push(`\\mathrm{${P_or_Sq}}^{${P}}`);
    }
    if(bocksteins >> Ps.length){
        simple_result.push("b");
        latex_result.push("\\beta");
    }
    simple_result = simple_result.join(" ");
    latex_result = latex_result.join(" ");
    return {
        "simple" : simple_result,
        "latex"  : latex_result
    };
}


function ademElementToString(algebra, degree, result_vector){
    let length = cAdemAlgebra_element_toString(algebra.string_buffer, algebra.c_algebra, degree, result_vector);
    let string_array = Module.HEAPU8.subarray(algebra.string_buffer, algebra.string_buffer + length);
    let result = Utf8ArrayToStr(string_array);
    let Sq_or_P = algebra.p === 2 ? "Sq" : "P";
    let latex_result = result.replace(/P(\d*)/g, `\\mathrm{${Sq_or_P}}^{$1}`);
    let simple_result = result.replace(/P/g,"");
    // let latex_result = [];
    // let simple_result = [];
    // for(let i=0; i<result_array.length; i++){
    //     let coeff = result_array[i]; 
    //     if(coeff === 0){
    //         continue;
    //     }
    //     let c_basis_element = cAdemAlgebra_basisElement_fromIndex(algebra.c_algebra, degree, i);
    //     let bocksteins = cAdemAlgebra_basisElement_getBocksteins(c_basis_element);
    //     let P_length = cAdemAlgebra_basisElement_getPlength(c_basis_element);
    //     let c_Ps_pointer = cAdemAlgebra_basisElement_getPs(c_basis_element);
    //     let c_Ps_array = Module.HEAPU32.subarray(c_Ps_pointer/sizeof_uint, c_Ps_pointer/sizeof_uint + P_length);
    //     let output = jsAdemMonomialToOutput(c_Ps_array, bocksteins, algebra.p);
    //     let latex_output = output.latex;
    //     let simple_output = output.simple;
    //     if(coeff !== 1){
    //         latex_output = `${coeff}*${latex_output}`;
    //         simple_output = `${coeff}*${simple_output}`;
    //     }
    //     latex_result.push(latex_output);
    //     simple_result.push(simple_output);
    // }
    // latex_result = latex_result.join(" + ");
    // simple_result = simple_result.join(" + ");
    // if(latex_result.length === 0){
    //     latex_result = "0";
    // }
    return {"latex" : latex_result, "simple" : simple_result};
}

handlers["construct"] = function construct(data){
    algebra = constructAlgebra(data.algebra);
}

handlers["calculate_adem"] = function compute(data){
    if(!algebra){
        Module.print("Assertion error: algebra isn't defined. I need to be called with construct first.");
        return;
    }
    if(!algebra.type == "adem"){
        Module.print("My algebra is not of type 'adem' =/");
        return;
    }
    let input = data.input;
    let p = algebra.p;
    let q = algebra.generic ? 2*p - 2 : 1;
    let jsMonomial = ademSequenceTojsAdemMonomial(q, input);
    if(jsMonomial.error){
        self.postMessage({
            "cmd" : "adem_parse_error",
            "error_str" : jsMonomial.error,
            "position" : jsMonomial.position
        });
        return;
    }
    let Ps = jsMonomial.Ps;
    let bocksteins = jsMonomial.bocksteins;
    let latex_input = jsAdemMonomialToOutput(jsMonomial.Ps, jsMonomial.bocksteins, algebra.p).latex;
    let degree = jsMonomial.degree;
    if(degree >= degree_caps[p]){
        let output_degree = degree;
        let max_degree = degree_caps[p];
        if(p==2 && bocksteins == 0){
            output_degree = Math.floor(output_degree/2);
            max_degree = Math.floor(max_degree/2);
            q = 1;
        }

        self.postMessage({
            "cmd" : "adem_error",
            "error_str" : `Degree ${output_degree} of input is too large. We run out of memory past degree ${max_degree}.`
        });
        return;
    }
    generateBasis(algebra, degree + 1);
    let c_Ps = Module._malloc(sizeof_uint * Ps.length);
    Module.HEAPU32.set(new Uint32Array(Ps), c_Ps/sizeof_uint);
    Module.print("degree:",degree)
    let c_basis_element = cAdemAlgebra_basisElement_construct(degree, Ps.length, c_Ps, bocksteins);
    let coeff = 1;
    let excess = 1<<20;
    let dimension = cAdemAlgebra_getDimension(algebra.c_algebra, degree, excess);
    let c_result_vector = cVector_construct(algebra.p, dimension, 0);
    cAdemAlgebra_makeMonoAdmissible(algebra.c_algebra, c_result_vector, coeff, c_basis_element, excess);
    // let c_result_pointer = Module._malloc(sizeof_uint * dimension);
    // cVector_unpack(c_result_pointer, c_result_vector);
    // let result_array = Module.HEAPU32.subarray(c_result_pointer/sizeof_uint, c_result_pointer/sizeof_uint + dimension);
    let result = ademElementToString(algebra, degree, c_result_vector);
    cVector_free(c_result_vector);
    self.postMessage({
        "cmd" : "adem_result", 
        "latex_input" : latex_input,
        "latex_result" : result["latex"], 
        "simple_result" : result["simple"]
    })
}


// https://stackoverflow.com/a/22373197/1942152
function Utf8ArrayToStr(array) {
    var out, i, len, c;
    var char2, char3;

    out = "";
    len = array.length;
    i = 0;
    while(i < len) {
    c = array[i++];
    switch(c >> 4)
    { 
      case 0: case 1: case 2: case 3: case 4: case 5: case 6: case 7:
        // 0xxxxxxx
        out += String.fromCharCode(c);
        break;
      case 12: case 13:
        // 110x xxxx   10xx xxxx
        char2 = array[i++];
        out += String.fromCharCode(((c & 0x1F) << 6) | (char2 & 0x3F));
        break;
      case 14:
        // 1110 xxxx  10xx xxxx  10xx xxxx
        char2 = array[i++];
        char3 = array[i++];
        out += String.fromCharCode(((c & 0x0F) << 12) |
                       ((char2 & 0x3F) << 6) |
                       ((char3 & 0x3F) << 0));
        break;
    }
    }

    return out;
}