from steenrod_module import *

module_list = []

S2 = FiniteSteenrodModule(p=2)
S2.name = "$S_2$"
S2.file_name = "S_2"
S2.x0 = S2.add_basis_element("x0", 0)
module_list.append(S2)

C2 = FiniteSteenrodModule(p=2)
C2.name = "$C(2)$"
C2.file_name = "C2"
C2.x0 = C2.add_basis_element("x0", 0)
C2.x1 = C2.add_basis_element("x1", 1)
C2.add_Sq_action(1, C2.x0, C2.x1)
C2.validate()
module_list.append(C2)

Ceta = FiniteSteenrodModule(p=2)
Ceta.name = "$C(\\eta)$"
Ceta.file_name = "Ceta"
Ceta.x0 = Ceta.add_basis_element("x0", 0)
Ceta.x2 = Ceta.add_basis_element("x2", 2)
Ceta.add_Sq_action(2, Ceta.x0, Ceta.x2)
Ceta.validate()
module_list.append(Ceta)

Cnu = FiniteSteenrodModule(p=2)
Cnu.name = "$C(\\nu)$"
Cnu.file_name = "Cnu"
Cnu.x0 = Cnu.add_basis_element("x0", 0)
Cnu.x4 = Cnu.add_basis_element("x4", 4)
Cnu.add_Sq_action(4, Cnu.x0, Cnu.x4)
Cnu.validate()
module_list.append(Cnu)

Csigma = FiniteSteenrodModule(p=2)
Csigma.name = "$C(\\sigma)$"
Csigma.file_name = "Csigma"
Csigma.x0 = Csigma.add_basis_element("x0", 0)
Csigma.x8 = Csigma.add_basis_element("x8", 8)
Csigma.add_Sq_action(8, Csigma.x0, Csigma.x8)
Csigma.validate()
module_list.append(Csigma)

C2_eta = FiniteSteenrodModule(p=2)
C2_eta.name = "$C(2,\\eta)$"
C2_eta.file_name = "C2_eta"
C2_eta.x0 = C2_eta.add_basis_element("x0", 0)
C2_eta.x1 = C2_eta.add_basis_element("x1", 1)
C2_eta.x3 = C2_eta.add_basis_element("x3", 3)
C2_eta.add_Sq_action(1, C2_eta.x0, C2_eta.x1)
C2_eta.add_Sq_action(2, C2_eta.x1, C2_eta.x3)
C2_eta.validate()
module_list.append(C2_eta)

S3 = FiniteSteenrodModule(p=3)
S3.name = "$S_3$"
S3.file_name = "S_3"
S3.x0 = S3.add_basis_element("x0", 0)
module_list.append(S3)

C3 = FiniteSteenrodModule(p=3)
C3.name = "$C(3)$"
C3.file_name = "C3"
C3.x0 = C3.add_basis_element("x0", 0)
C3.x1 = C3.add_basis_element("x1", 1)
C3.add_b_action(C3.x0, C3.x1)
C3.validate()
module_list.append(C3)

Calpha = FiniteSteenrodModule(p=3)
Calpha.name = "$C(\\alpha)$"
Calpha.file_name = "Calpha"
Calpha.x0 = Calpha.add_basis_element("x0", 0)
Calpha.x4 = Calpha.add_basis_element("x4", 4)
Calpha.add_P_action(1, Calpha.x0, Calpha.x4)
Calpha.validate()
module_list.append(Calpha)


Calphaalpha = FiniteSteenrodModule(p=3)
Calphaalpha.name = "$C(\\alpha, \\alpha)$"
Calphaalpha.file_name = "X3"
Calphaalpha.x0 = Calphaalpha.add_basis_element("x0", 0)
Calphaalpha.x4 = Calphaalpha.add_basis_element("x4", 4)
Calphaalpha.x8 = Calphaalpha.add_basis_element("x8", 8)
Calphaalpha.add_P_action(1, Calphaalpha.x0, Calphaalpha.x4)
Calphaalpha.add_P_action(1, Calphaalpha.x4, Calphaalpha.x8)
Calphaalpha.add_P_action(2, Calphaalpha.x0, 2*Calphaalpha.x8)
Calphaalpha.validate()
module_list.append(Calphaalpha)


Joker = FiniteSteenrodModule(p=2)
Joker.name = "Joker"
Joker.file_name = "Joker"
x0 = Joker.add_basis_element("x0", 0)
x1 = Joker.add_basis_element("x1", 1)
x2 = Joker.add_basis_element("x2", 2)
x3 = Joker.add_basis_element("x3", 3)
x4 = Joker.add_basis_element("x4", 4)
Joker.add_Sq_action(2, x0, x2)
Joker.add_Sq_action(2, "x2", x4)
Joker.add_Sq_action(1, "x0", x1)
Joker.add_Sq_action(2, "x1", x3)
Joker.add_Sq_action(1, "x3", x4)
Joker.add_Sq_action(3, "x1", x4)
Joker.validate()
module_list.append(Joker)

RP4 = FiniteSteenrodModule.from_Bruner_file("modules/RP4")
RP4.name = "$\\mathbb{RP}^4$"
RP4.file_name = "RP4"
RP4.validate()
module_list.append(RP4)

ko = FiniteSteenrodModule(p=2)
ko.name = "$\\mathit{ko}$"
ko.file_name = "ko"
ko.add_basis_element("x0", 0)
ko.profile = {"truncated" : True, "p_part" : [2,1]}
module_list.append(ko)

tmf2 = FiniteSteenrodModule(p=2)
tmf2.name = "$\\mathit{tmf}_{(2)}"
tmf2.file_name = "tmf2"
tmf2.add_basis_element("x0", 0)
tmf2.profile = {"truncated" : True,"p_part" : [3,2,1]}
module_list.append(tmf2)


def standard_modules_to_json():
    for M in module_list:
        print("modules/%s.json" % M.file_name)
        M.to_json_file("modules/%s.json" % M.file_name)
    
    module_registry = [
        {
            "p" : M.p, 
            "generic" : M.generic, 
            "name" : M.name, 
            "file_name" : M.file_name
        } 
        for M in module_list
    ]
    json_str = json.dumps(module_registry)
    write_file("modules/module_registry.json", json_str)

if __name__ == "__main__":
    standard_modules_to_json()