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
        return super().read(cls.path + cls.file)
    
    def save_update(self):
        
        self._update_log()
        self.write(self.path + self.file, overwrite=True)
        
    def save_to_html(self):
        
        self.write(self.path + self.vfile, format='jsviewer')
        
    def _update_log(self):
        
        with open(self.path+self._logfile, 'a+') as f:

            f.write('----------------------\n')
            f.write(str(datetime.datetime.now())+'\n')
            f.write('MODIFICATION BY {} :\n'.format(getuser()))
            f.writelines('\n'.join(self.log)+'\n')
            f.close()
    
#     def __init__(self, target_list=None, verbose=True):
        
#         if verbose: print('Reading masterfile.csv...')
#         all_data = Table().read(self.path+self.file, format='ascii.ecsv')
#         all_data.main_col = self.main_col

#         if target_list:
#             super().__init__(all_data.by_plName(*target_list))
#         else:
#             super().__init__(all_data)