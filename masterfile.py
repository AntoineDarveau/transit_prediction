from transit_prediction.utils.table_edited import Table
from getpass import getuser
import datetime

class MasterFile(Table):
    
    path = "/home/adb/Archive/"
    file = "masterfile.ecsv"
    vfile = "masterfile.html"
    main_col = 'pl_name'
    _logfile = 'log_masterfile_temp.txt'
    
    @classmethod
    def read(cls):
        out = super().read(cls.path + cls.file)
        cls.save_dim = (len(out), len(out.keys()))
        return out
    
    def save_update(self):
        new_dim = len(self), len(self.keys())
        if (self.save_dim[0]<=new_dim[0] 
         and self.save_dim[1]<=new_dim[1]):
            self._update_log()
            self.write(self.path + self.file, overwrite=True)
        else:
            raise ValueError("dimensions don't match \n"
                            +"masterfile dim: {} \n".format(self.save_dim)
                            +"table dim: {}".format(new_dim))
    def save_to_html(self):
        
        self.write(self.path + self.vfile, format='jsviewer')
        
    def _update_log(self):
        
        with open(self.path+self._logfile, 'a+') as f:

            f.write('----------------------\n')
            f.write(str(datetime.datetime.now())+'\n')
            f.write('MODIFICATION BY {} :\n'.format(getuser()))
            f.writelines('\n'.join(self.log)+'\n')
            f.close()
        self.log = []  # Erase log
    
#     def __init__(self, target_list=None, verbose=True):
        
#         if verbose: print('Reading masterfile.csv...')
#         all_data = Table().read(self.path+self.file, format='ascii.ecsv')
#         all_data.main_col = self.main_col

#         if target_list:
#             super().__init__(all_data.by_plName(*target_list))
#         else:
#             super().__init__(all_data)