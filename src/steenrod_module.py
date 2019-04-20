import json
import itertools
import math
from functools import reduce

from infinity import Infinity
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
        super(SteenrodModuleElement, self).__init__(module.p, dictionary)
        
    def adem_act(self, adem_elt):
        steenrod.implementedByAssignmentLaterInThisFile()

    def milnor_act(self, adem_elt):
        steenrod.implementedByAssignmentLaterInThisFile()
        
    def basis_degree(self, b):
        return self.module.gens[b]

class SteenrodModule:
    def __init__(self, *, p, name = None, generic = None, gens = None, sq_actions = None, milnor_actions = None):
        self.p = p
        self.name = name
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
        self.validated = True
        self.failed_relations = {} 
        self.milnor_actions = sq_actions or []
    
    def __repr__(self):
        args = []
        args += ["p = %s" % self.p]
        if self.p == 2 and self.generic:
            args += ["generic"]
        args += [ "%s generators" % len(self.gens) ]
        if not self.validated:
            args += ["not validated"]
        elif len(self.failed_relations) > 0:
            num = len(self.failed_relations)
            fails_num_str = "fails %s relation%s" % (num, "s" if num > 1 else "")
            args += [fails_num_str]
        args = ", ".join(args)
        result = ""
        if self.name:
            result += self.name + " : "
        result += "SteenrodModule(%s)" % args
        return result
    
    #
    # Some elements 
    #
    
    def getBasisElement(self, b):
        return SteenrodModuleElement({b : 1}, module = self)
        
    def getElement(self, dict):
        return SteenrodModuleElement(dict, module = self)
        
    def zero(self):
        return SteenrodModuleElement({}, module = self)
    
    #
    # Computing the action of a general Adem element from lookup table
    #
    
    def sq_act_on_basis(self, module_basis_elt, sq):
        """
            Used for adem_basis_act which is extended to the general action
        """
        if sq == 0:
            return self.module.getBasisElement( module_basis_elt )
        actions = self.sq_actions
        key = (sq, module_basis_elt)
        return actions[key] if key in actions else {}
        
    def sq_act_on_vector(self, module_element, sq):
        """
            Used for adem_basis_act which is extended to the general action.
        """
        if sq == 0:
            return module_vector
        actions = self.sq_actions
        result = module_element.module.zero()
        for (module_basis_elt, coeff) in module_element.items():
            output = self.sq_act_on_basis(module_basis_elt, sq)
            result.add_in_place(output, coeff)
        return result
                 
    @staticmethod
    def adem_basis_act(module_basis_elt, adem_basis_elt, *, module):
        """
            This is for handing to linearly_extend_map. The argument module will be filled by 
            basis_elt.algebra which will be a Steenrod_Module. 
            algebra -- will be this.
        """
        sqs = adem_basis_elt
        if module.generic:
            sqs = adem.adem_basis_elt_generic_map(P_fn = lambda P : P, b = 'b', basis_elt = adem_basis_elt)
        sqs = sqs[::-1] # Reverse
        module_vector = module.getBasisElement(module_basis_elt)
        result = reduce(module.sq_act_on_vector, sqs, module_vector)
        return result
        
    def generate_milnor_action(self):
        raise NotImplementedError()
        if not self.validQ():
            raise ValueError("Validate module first.")
        
    
    @staticmethod
    def milnor_basis_act(basis_elt, adem_elt, *, module):
        """
            This is for handing to linearly_extend_map. The argument module will be filled by 
            basis_elt.module which will be a Steenrod_Module. 
            module -- will be this.
        """
        raise NotImplementedError()
        #return module.milnor_actions[]
    
    
    #
    # Adding generators and relations
    #
        
    def add_basis_element(self, name, degree):
        self.gens[name] = degree
        elt = self.getBasisElement(name)
        self.sq_actions[(0, name)] = elt
        return elt
        
    def add_Sq_action(self, Sq, input_gen, output_vec):
        input_deg = self.gens[input_gen]
        self.check_degree(Sq, input_gen, output_vec)
        
        if self.validated:
            self.clear_validation()
        key = (Sq, input_gen)
        if key in self.sq_actions:
            output_vec += self.sq_actions[key]
        output_vec = self.getElement(output_vec)            
        self.sq_actions[key] = output_vec
        
    def add_P_action(self, P, input_gen, output_vec):
        if P == 'b':
            self.add_b_action(input_gen, output_vec)
            return 
        self.check_degree(2*(self.p-1)*P, input_gen, output_vec)
        if self.validated:
            self.clear_validation()
        key = (P, input_gen)
        if key in self.sq_actions:
            output_vec += self.sq_actions[key]
        output_vec = self.getElement(output_vec)             
        self.sq_actions[key] = output_vec
    
    def add_b_action(self, input_gen, output_vec):
        self.check_degree(1, input_gen, output_vec)
        if self.validated:
            self.clear_validation()
        key = ('b', input_gen)
        if key in self.sq_actions:
            output_vec += self.sq_actions[key]
        output_vec = self.getElement(output_vec)             
        self.sq_actions[key] = output_vec
        
    def add_milnor_Q_action(self, Q, input_gen, output_vec):
        raise NotImplementedError()
    
    def add_milnor_P_action(self, Q, input_gen, output_vec):
        raise NotImplementedError()
        
    def check_degree(self, operator_degree, input_gen, output_vec):
        input_degree  = self.gens[input_gen]
        output_degree = output_vec.degree()
        if output_degree is None:
            raise ValueError("Output vector %s is not homogenous" % output_vec)
        if operator_degree + input_degree != output_degree:
            raise ValueError("Invalid relation: degree of output vector should be %s but instead is %s" 
                % (operator_degree + input_degree, output_degree))
    
    def sq_action_strs(self):
        result = []
        for ((sq, input), output) in [ x for x in self.sq_actions.items() if x[0][0] != 0 ]:
            if generic:
                op = "b" if sq == "b" else "P" + sq
            else:
                op = "Sq" + sq
            result += ["%s(%s) = %s" % (op, input, output)]
        return result
        
    def clear_validation(self):
        self.validated = False
        self.failed_relations = {}
    
    def validate_generic(self):
        self.failed_relations = {}
        self.validated = True
        max_deg = max(self.gens.values())
        p = self.p
        q = 2*(p - 1) 
        for gen, gen_deg in self.gens.items():
            for relation_dim in range(2, (max_deg - gen_deg)//q + 1):
                for epsilon in [0, 1]:
                    # We want Pi*b*Pj inadmissible so that means i < p * j + epsilon.
                    # relation_dim = i + j so j = relation_dim - i
                    # so i < p * (relation_dim - i) + epsilon so i < (p * relation_dim + epsilon) / ( p + 1)
                    # We need to round up so that Python includes the last integer 
                    # if 2 * relation_dim / 3 is not an integer
                    for i in range(1, int(math.ceil((p * relation_dim + epsilon) / ( p + 1)))):
                        j = relation_dim - i
                        Pi = self.adem_algebra.P(i)
                        Pj = self.adem_algebra.P(j)
                        b = self.adem_algebra.b_or_unit(epsilon)
                        v = self.getBasisElement( gen )
                        boundary = (Pi * b * Pj) * v - Pi * (b * Pj * v)
                        if str(boundary) != '0':
                            self.failed_relations[(i, epsilon, j, gen)] = boundary
        return self
    
    def validate_2(self):
        self.failed_relations = {}
        self.validated = True
        max_dim = max(self.gens.values())
        for gen, degree in self.gens.items():
            for relation_dim in range(2, max_dim - degree + 1):
                # We want Sqi*Sqj inadmissible so that means i < 2 * j.
                # relation_dim = i + j so j = relation_dim - i
                # so i < 2 * (relation_dim - i) so i < 2 * relation_dim / 3. 
                # We need to round up so that Python includes the last integer 
                # if 2 * relation_dim / 3 is not an integer
                for i in range(1, int(math.ceil((2 * relation_dim) / 3))):
                    j = relation_dim - i
                    Sqi = self.adem_algebra.Sq(i)
                    Sqj = self.adem_algebra.Sq(j)
                    v = self.getBasisElement( gen )
                    boundary = (Sqi * Sqj) * v - Sqi * (Sqj * v)
                    if str(boundary) != '0':
                        self.failed_relations[(i, j, gen)] = boundary
        return self
        
    def validate(self):
        if self.generic:
            self.validate_generic()
        else:
            self.validate_2()
        return self
                        
    def getFailedRelationStrings_generic(self):
        relation_strings = []
        for reln, boundary in self.failed_relations.items():
            i = reln[0]
            epsilon = reln[1]
            j = reln[2]
            v = reln[3]
            Pi = self.adem_algebra.P(i)
            Pj = self.adem_algebra.P(j)
            b_or_unit = self.adem_algebra.b_or_unit(epsilon)
            relation = Pi * b_or_unit * Pj
            Pj_or_bPj_str  = "b " if epsilon else ""
            Pj_or_bPj_str += str(Pj)
            rel_str = "(%s) * %s  -  %s * (%s * %s) = %s" % (relation, v,      Pi, Pj_or_bPj_str, v,  boundary)
            relation_strings += [rel_str]
        return relation_strings                        
                        
    def getFailedRelationStrings_2(self):
        relation_strings = []
        for reln, boundary in self.failed_relations.items():
            i = reln[0]
            j = reln[1]
            v = reln[2]
            Sqi = self.adem_algebra.Sq(i)
            Sqj = self.adem_algebra.Sq(j)
            rel_str = "(%s) * %s  -  %s * (%s * %s) = %s" % (Sqi * Sqj, v, Sqi, Sqj, v, boundary)
            relation_strings += [rel_str]
        return relation_strings
    
    def getFailedRelationStrings(self):
        if self.generic:
            return self.getFailedRelationStrings_generic()
        else:
            return self.getFailedRelationStrings_2()
    
    def validQ(self):
        return self.validated and len(self.failed_relations) == 0
     
    #
    # Making new modules
    #   

    ##### @staticmethod
    def tensor(*modules):
        """
            Can be used either as a static method "SteenrodModule.tensor(M, N, P)" or like "M.tensor(N)".
            For "M.tensor(N)" we offer the syntactic sugar M*N. Pretty sweet =)
        """
        p = modules[0].p
        generic = modules[0].generic
        for M in modules:
            if not M.validQ():
                raise ValueError("We only tensor validated modules. Run M.validate() first.")
            if M.p != p:
                raise ValueError("We only tensor modules that share the same prime.")
            if M.generic != generic:
                raise ValueError("We only tensor modules that share the same genericness.")
            
        result = SteenrodModule(p = p, generic = generic)
        for x in itertools.product(*[M.gens.items() for M in modules]):
            (gens, dims) = zip(*x)
            tensor = Vector.tensor_basis_elements(*gens)
            dim = sum(dims)
            result.add_basis_element(tensor, dim)
        
        if generic:
            add_action = result.add_P_action
        else:
            add_action = result.add_Sq_action
        
        for sq_input_output in itertools.product(*[[kv for kv in M.sq_actions.items() if kv[0][0] != 'b'] for M in modules]):
            # Zip is transpose. We want the sqs, inputs, and outputs each in their own list
            (sq_input, outputs) = zip(*sq_input_output)
            (sqs, inputs) = zip(*sq_input)
            sq = sum(sqs)
            if sq == 0:
                continue
            input = Vector.tensor_basis_elements(*inputs)
            output = Vector.tensor(*outputs, output_module = result )
            add_action(sq, input, output)
            
        if generic:
            for i in range(0,len(modules)):
                pre = modules[:i]
                M = modules[i]
                post = modules[i+1:]
                pre_elts  = [[(b,result.getBasisElement(b)) for b in N.gens] for N in pre]
                post_elts = [[(b,result.getBasisElement(b)) for b in N.gens] for N in post]
                act = [(input,v) for ((sq, input),v) in M.sq_actions.items() if sq == 'b']
                for input_output in itertools.product(*pre_elts, act, *post_elts):
                    (inputs, outputs) = zip(*input_output)
                    input = Vector.tensor_basis_elements(*inputs)
                    output = Vector.tensor(*outputs, output_module = result )
                    result.add_b_action(input, output)
        
        result.validated = True
        return result
            
    def dualize(self):
        raise NotImplementedError()
        result = SteenrodModule(p = p, generic = generic)
        return result        
        
    def truncate(self, min_dim = -Infinity, max_dim = Infinity):
        if min_dim is None:
            min_dim = -Infinity
        if max_dim is None:
            max_dim = Infinity

        p = self.p
        generic = self.generic
        result = SteenrodModule(p = p, generic = generic)
        for (gen, deg) in self.gens.items():
            if(min_dim <= deg <= max_dim):
                result.add_basis_element(gen, deg)
        
        for ((sq, input), output) in self.sq_actions.items():
            op_dim = sq
            if generic and sq == 'b':
                op_dim = 1
            elif generic:
                op_dim *= 2*(self.p - 1)
            add_action = result.add_P_action if generic else result.add_Sq_action
            input_dim = self.gens[input]
            output_dim = input_dim + op_dim
            if(min_dim <= input_dim and output_dim <= max_dim):
                add_action(sq, input, output)
        result.validated = True
        return result
            
            
    # Some syntactic sugar  
    def __mul__(self, other):
        """
            M*N is the tensor product of M and N
        """
        return SteenrodModule.tensor(self,other)
    
    def __invert__(self):
        """
            ~M is M dual
        """
        return self.dualize()
        
    def __getitem__(self, key):
        """
            Take all dimension between min and max inclusive.
            Not consistent with normal Python slicing which is inclusive of min
            and exclusive of max, but I don't care.
        """
        if isinstance( key, slice ) :
            return self.truncate(key.start, key.stop)
        else:
            raise ValueError("You must give a slice M[min::max]")

    #
    # Serialization and deserialization.
    #
    
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
    A = M.adem_algebra
    Sq1 = A.Sq(1)
    Sq2 = A.Sq(2)
    Sq3 = A.Sq(3)
    x0 = M.add_basis_element("x0",0)
    x1 = M.add_basis_element("x1",1)
    x2 = M.add_basis_element("x2",2)
    x3 = M.add_basis_element("x3",3)
    x4 = M.add_basis_element("x4",4)
    M.add_Sq_action(2, "x0", x2)
    M.add_Sq_action(2, "x2", x4)
    M.add_Sq_action(1, "x0", x1)
    M.add_Sq_action(2, "x1", x3)
    M.add_Sq_action(1, "x3", x4)
    M.add_Sq_action(3, "x1", x4)  
    M.validate()
    M3 = SteenrodModule(p = 3)
    A3 = M3.adem_algebra
    P1 = A3.P(1)
    P2 = A3.P(2)
    P3 = A3.P(3)
    y0 = M3.add_basis_element("y0",0)
    y1 = M3.add_basis_element("y1",1)
    y4 = M3.add_basis_element("y4",4)
    y5 = M3.add_basis_element("y5",5)
    M3.add_P_action(1, "y0", y4)
    M3.add_P_action(1, "y1", y5)
    M3.add_b_action("y0", y1)
    M3.add_b_action("y4", y5)
    M3.validate()
    
    
