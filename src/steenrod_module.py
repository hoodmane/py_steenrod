import json
import itertools
import math
from functools import reduce

from FpVectorSpace import Vector
import adem
import steenrod

def read_file(input_file):
    f = open(input_file, 'r')
    contents = f.read()
    f.close()
    return contents
    
def write_file(output_file, contents):
    f = open(output_file, 'w')
    f.write(contents)
    f.close()    

class SteenrodModuleElement(Vector):
    def __init__(self, dictionary = None, *, module):
        self.module = module
        self.algebra = module
        super(SteenrodModuleElement, self).__init__(module.p, dictionary)
        
    def adem_act(self, adem_elt):
        steenrod.implementedByAssignmentLaterInThisFile()

    def milnor_act(self, adem_elt):
        steenrod.implementedByAssignmentLaterInThisFile()    

class SteenrodModule:
    def __init__(self, *, p, generic = None, gens = None, sq_actions = None, milnor_actions = None):
        self.p = p
        if generic is None:
            generic = p != 2
        self.generic = generic
        self.adem_algebra = steenrod.AdemAlgebra.getInstance( self.p, self.generic )
        self.gens = gens or {}
        sq_actions = sq_actions or []
        self.sq_actions = {}
        for act in sq_actions:
            if "Sq" in act:
                self.add_Sq_action(self, act["Sq"], act["input_gen"], act["output_vec"])
            elif "P" in act:
                self.add_P_action(self, act["P"], act["input_gen"], act["output_vec"])
            elif "b" in act:
                self.add_b_action(self, act["input_gen"], act["output_vec"])
            else:
                raise ValueError()
        self.verified = True
        self.failed_relations = {} 
        self.milnor_actions = sq_actions or []
    
    def getBasisElement(self, b):
        return SteenrodModuleElement({b : 1}, module = self)
        
    def zero(self):
        return SteenrodModuleElement({}, module = self)
    
    def sq_act(self, module_basis_elt, sq):
        """
            Used for adem_basis_act which is extended to the general action
        """
        if sq == 0:
            return algebra.getBasisElement( module_basis_elt )
        actions = self.sq_actions
        key = (sq, module_basis_elt)
        return actions[key] if key in actions else None
                 
    @staticmethod
    def adem_basis_act(module_basis_elt, adem_basis_elt, *, algebra):
        """
            This is for handing to linearly_extend_map. The argument algebra will be filled by 
            basis_elt.algebra which will be a Steenrod_Module. Maybe this argument is misnamed...
            algebra -- will be this.
        """
        sqs = adem_basis_elt
        if algebra.generic:
            sqs = adem.adem_basis_elt_generic_map(P_fn = lambda P : P, b = 'b', basis_elt = adem_basis_elt)
        result = reduce(algebra.sq_act, sqs, module_basis_elt)
        #print(adem_basis_elt, module_basis_elt, result)
        return result
        
    @staticmethod
    def milnor_basis_act(basis_elt, adem_elt, *, algebra):
        """
            This is for handing to linearly_extend_map. The argument algebra will be filled by 
            basis_elt.algebra which will be a Steenrod_Module. Maybe this argument is misnamed...
            algebra -- will be this.
        """
        raise NotImplementedError()
        #return algebra.milnor_actions[
    
        
    def add_generator(self, name, degree):
        self.gens[name] = degree
        
    def add_Sq_action(self, Sq, input_gen, output_vec):
        self.clear_verification()
        self.sq_actions[(Sq, input_gen)] = output_vec
        
    def add_P_action(self, P, input_gen, output_vec):
        self.clear_verification()
        self.sq_actions[(P, input_gen)] = output_vec
    
    def add_b_action(self, input_gen, output_vec):
        self.clear_verification()
        self.sq_actions[('b', input_gen)] = output_vec
        
    def add_milnor_Q_action(self, Q, input_gen, output_vec):
        raise NotImplementedError()
    
    def add_milnor_P_action(self, Q, input_gen, output_vec):
        raise NotImplementedError()
        
    def clear_verification(self):
        self.verified = False
        self.failed_relations = {}
    
    def verify(self):
        self.clear_verification()
        
        max_dim = max(self.gens.values())
        for gen, degree in self.gens.items():
            for dimop in range(2, max_dim - degree + 1):
                for a in range(1, int(math.ceil((2 * dimop) / 3))):
                    b = dimop - a
                    Sqa = self.adem_algebra.Sq(a)
                    Sqb = self.adem_algebra.Sq(b)
                    v = self.getBasisElement( gen )
                    boundary = (Sqa * Sqb) * v - Sqa * (Sqb * v)
                    if str(boundary) != '0':
                        self.failed_relations[(a, b, gen)] = boundary
                        
    def getFailedRelationStrings(self):
        relation_strings = []
        for reln, boundary in self.failed_relations.items():
            a = reln[0]
            b = reln[1]
            v = reln[2]
            Sqa = self.adem_algebra.Sq(2)
            Sqb = self.adem_algebra.Sq(2)
            rel_str = "(%s) * %s  -  %s * (%s * %s) = %s" % (Sqa*Sqb, v, Sqa, Sqb, v, boundary)
            relation_strings += [rel_str]
        return relation_strings

    @staticmethod
    def tensor_product(M1, M2):
        raise NotImplementedError()
    
            
    def dualize(self):
        raise NotImplementedError()
    
    def __mul__(self, other):
        return SteenrodModule.tensor_product(self,other)
    
    # This overloads ~
    def __invert__(self):
        return self.dualize()
        
    @staticmethod
    def from_Bruner_str(str):
        lines = [s for s in module.split("\n") if s]
        num_gens = int(lines[0])
        gens =  [int(s) for s in lines[1].split(" ") if s]
        gens_dict = {}
        idx = 0
        for g in gens:
            gens_dict[str(idx)] = g
            idx += 1
        actions = []
        for l in lines[2:]:
            l = [int(s) for s in l.split(" ") if s]
            act = {"Sq" : l[1], "input" : str(l[0]), "output" : []}
            for e in l[3:]:
                act["output"].append([1, str(e)])
            actions.append(act)
        return {"generators" : gens_dict, "actions" : actions }
        
    def to_Bruner_str(self):
        gens = self.generators
        actions = self.sq_actions
        num_gens = str(len(gens))
        gen_degrees = " ".join([str(x) for x in gens.values()])
        # We need a map (index in generator list) --> (generator name) because Bruner modules only refer to elts by their index.
        l = [(gens[k], k) for k in gens]
        l.sort()
        gen_name_to_idx = { str(ind) : x[1] for ind, x in enumerate(l)}
        output_lines = [num_gens]
        output_lines += [gen_degrees]
        for act in actions:
            input_idx = gen_name_to_idx[act["input"]]
            sq = str(act["Sq"])
            output_length = str(len(act["output"]))
            if generic:
                output_list = act["output"]
                coeffs        = [str(x[0]) for x in output_list]
                basis_indices = [str(gen_name_to_idx[x[1]]) for x in output_list]
                output_vector = [None] * (2*len(output_list))
                output_vector[::2] = coeffs
                output_vector[1::2] = basis_indices 
            else:
                output_vector = [ str(gen_name_to_idx[x[1]]) for x in act["output"] ]
            action_list   = [input_idx, sq, output_length ] 
            action_list  += output_vector            
            output_lines += [ " ".join(action_list) ]
        return "\n".join(output_lines)
        
    @staticmethod
    def from_json_obj(obj):
        raise NotImplementedError()
        
    def to_json_obj(self):
        raise NotImplementedError()        

    @staticmethod
    def from_Bruner_file(file):    
        file_contents = read_file(file)
        return SteenrodModule.from_Bruner_str(file_contents)

    def to_Bruner_file(self, file):
        write_file(file, self.to_Bruner_str())

    @staticmethod
    def from_json_file(file):    
        json = json.loads(read_file(file))
        return SteenrodModule.from_json_obj(json)  

    def to_json_file(self, file):
        obj = self.to_json_obj()
        json_str = json.dumps(obj)
        write_file(file, json_str)


SteenrodModuleElement.adem_act = Vector.linearly_extend_map(SteenrodModule.adem_basis_act)


if __name__ == "__main__":
    M = SteenrodModule(p = 2)
    M.add_generator("x0",0)
    M.add_generator("x2",2)
    M.add_generator("x4",4)
    M.add_Sq_action(2, "x0", M.getBasisElement("x2"))
    M.add_Sq_action(2, "x2", M.getBasisElement("x4"))
    M.verify()
