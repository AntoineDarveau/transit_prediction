

class List_of_constraints(list):
    
    def show(self):
        out = {}
        for el in self:
            name = type(el).__name__
            infos = vars(el)
            out[name] = infos
        return out