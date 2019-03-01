from transit_prediction.utils.table_edited import Table

class MasterFile(Table):
    
    path = "/home/adb/Archive/"
    file = "masterfile.ecsv"
    vfile = "masterfile.html"
    main_col = 'pl_name'
    
    @classmethod
    def read(cls):
        return super().read(cls.path + cls.file)
    
    def save_update(self):
        
        self.write(self.path + self.file, overwrite=True)
        self.write(self.path + self.vfile, format='jsviewer')
    
#     def __init__(self, target_list=None, verbose=True):
        
#         if verbose: print('Reading masterfile.csv...')
#         all_data = Table().read(self.path+self.file, format='ascii.ecsv')
#         all_data.main_col = self.main_col

#         if target_list:
#             super().__init__(all_data.by_plName(*target_list))
#         else:
#             super().__init__(all_data)