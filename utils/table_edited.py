from astropy.table import join
import astropy.table as table
from astropy.table.operations import _join, _merge_table_meta
import numpy as np
from collections import OrderedDict
from astropy.units import Unit, UnitTypeError 
import warnings
from astropy.utils.exceptions import AstropyUserWarning


class MaskedColumn(table.MaskedColumn):

    def find(self, sub, start=0, end=None):
        if isinstance(self, (MaskedColumn)):
            str_array = np.array(self, dtype=str)
            index = np.core.defchararray.find(str_array, sub, start=start, end=end) != -1
            return np.where(index)[0], str_array[index]

        else:
            return NotImplemented


class Column(table.Column):

    def find(self, sub, start=0, end=None):
        if isinstance(self, (Column)):
            str_array = np.array(self, dtype=str)
            index = np.core.defchararray.find(str_array, sub, start=start, end=end) != -1
            return np.where(index)[0], str_array[index]

        else:
            return NotImplemented


class Table(table.Table):

    # Redefine class attributes (if not, the originals would be taken)
    Column = Column
    MaskedColumn = MaskedColumn

    # Set attributes
    main_col = None  # Default column used to order
    log = []  # Save output when using insert_value method

    # New methods
    def rename_columns(self, old, new):

        for ko, kn in zip(old, new):
            self.rename_column(ko, kn)

    def nan_to_mask(self):
        """
        Replace nan by masked array
        """
        if self.masked:
            for k in self.keys():
                if self[k].dtype == float:
                    self[k].mask = np.isnan(self[k])
        else:
            raise TypeError("Input must be a Masked Table." +
                            "\n \t Set its mask to True before calling" +
                            " (example: t = Table(t,masked=True)).")

    def by_plName(self, *plName, name_key=None, remove=False):
        """
        Return the complete line of a given planet name (plName)
        """
        position = self.get_index(*plName, name_key=name_key)

        out = self[position]

        if remove and len(position) == 1:
            print(str(*plName) + ' has been removed')
            self.remove_row(int(position[0]))

        return out
    
    def get_index(self, *plName, name_key=None):
        '''
        Return the lines index where plName are located for the column given by name_key
        name_key default is given by main_col attribute of the object
        '''
        name_key = name_key or self.main_col

        position = []
        for pl in plName:
            try:
                position.append(int(self[name_key].find(pl)[0]))
            except TypeError:
                print('Wrong name or incomplete:')
                print(self[name_key].find(pl))
        
        return position
    
    def set_main_col(self, colname=None, extension='_temp'):
        '''
        Set self.main_col and assign it to the first column.
        If colname is None, simply assign self.main_col to the
        first column.
        '''
        if self.main_col is None:
            self.main_col = colname
        elif colname is None:
            colname = self.main_col
        
        colname_temp = colname+extension
        self.rename_column(colname, colname_temp)
        self.add_column(self[colname_temp], name=colname, index=0)
        self.remove_column(colname_temp)
    
    def correct_units(self, badunits=['degrees', 'days', 'hours','jovMass'],
                     gunits=['degree', 'day', 'hour','jupiterMass'], verbose=True,
                     debug=False):
        '''
        Correct columns units for astropy units
        '''
        text_frame = "Column {} corrected for '{}' unit (previous was '{}')"

        for col in self.colnames:
            if debug:
                print(col, self[col].unit)

            # Search for bad units
            for bunit, gunit in zip(badunits, gunits):
                if self[col].unit == bunit:
                    self[col].unit = gunit
                    
                    # Message and log it
                    self.log.append(
                        text_frame.format(col, self[col].unit, bunit))
                    if verbose:
                        print(self.log[-1])

    def cols_2_qarr(self, *keys):
        '''
        Returns columns given in input as astropy q_arrays
        '''
        out = []
        for k in keys:
            try:
                out.append(np.ma.array(self[k].data) * self[k].unit)
            except TypeError:
                out.append(np.ma.array(self[k].data))

        return tuple(out)

    def set_units(self, units, cols=None):
        '''
        Assign units to columns.
        units:
          list of units (str or astropy units) to be assign
        cols:
          list of columns names (str).
          If None is given, it takes all the keys, so Table.keys() as default
        '''

        if not cols:
            cols = self.keys()

        for col, u in zip(cols, units):
            self[col].unit = u

    def new_value(self, plName, col, value):

        names = np.array(self[self.main_col], dtype=str)
        position = np.where(names == plName)[0]

        self[col][position] = value

    def complete(self, right, key=None, join_type='left',
                 add_col=True, metadata_conflicts='warn', 
                 verbose=True, debug=False, **kwargs):
        """
        Add every missing data in self if present in right.

        join_type : 'inner': do not add new rows
                    'outer': add new rows if not present in self
        add_col: add new colums from right
        """

        key = key or self.main_col

        # Try converting inputs to Table as needed
        if not isinstance(right, Table):
            right = Table(right)

        # If not masked, simple join() will do
        if self.masked:
            out = self._complete(right, key=key, join_type=join_type,
                                 add_col=add_col, verbose=verbose, debug=debug)
        else:
            col_name_map = OrderedDict()
            out = _join(self, right, join_type, col_name_map, keys=key, **kwargs)

        # Merge the column and table meta data. Table subclasses might override
        # these methods for custom merge behavior.
        _merge_table_meta(out, [self, right], metadata_conflicts=metadata_conflicts)

        return out

    def _complete(self, right, key=None, join_type='left', add_col=True,
                  verbose=True, debug=False):

        if not key:
            raise ValueError('key is empty')

        # Save shared columns without "key"
        cols = intersection(self.keys(), right.keys())
        cols.remove(key)

        # Join tables
        join_t = join(self, right, join_type=join_type, keys=key)

        # Complete masked values of "self" if available in "right"
        for col in cols:

            # Add eventually a condition to check units!

            # Names of joined columns (default from join())
            col1, col2 = col + '_1', col + '_2'

            # Index of masked in "self" and not masked in "right"
            index = join_t[col1].mask & ~join_t[col2].mask

            # Reassign value
            join_t[col1].unshare_mask()
            join_t[col1][index] = join_t[col2][index]

            # Remove 2nd column and rename to original
            join_t[col1].name = col
            del join_t[col2]

        # Remove added columns from "right" if not wanted
        supp_cols = difference(right.keys(), self.keys())
        if debug: print(supp_cols)

        if not add_col and supp_cols:
            if verbose:
                print('remove non shared columns from second table')
            join_t.remove_columns(supp_cols)

        return join_t

    def complete_cols(self, col_in, col_out, name_key=None):
        '''
        Use a column from table to complete another column.
        Input:
            col_in: list of names of columns to use (list of str)
            col_out: list of names of columns to complete (list of str)
        '''
        # Take default col if none is given
        name_key = name_key or self.main_col

        # Def table with cols to use and rename it to cols to complete
        temp_table = Table(self[[name_key] + col_in], masked=True)
        temp_table.nan_to_mask()
        temp_table.rename_columns(col_in, col_out)

        # Complete with the temp_table
        return self.complete(temp_table, key=name_key)

    def add_calc_col(self, fct, *args, f_args=(), f_kwargs={}, col_keys=[], **kwargs):
        '''
        def add_calc_col(self, fct, *args, f_args=(), f_kwargs={}, col_keys=[], **kwargs):

        Add new column wich is the result of fct(table[col_keys], *f_args, **f_kwargs)

        args and kwargs are passed to MaskedColumn instantiation

        '''

        # Build tuple of columns inputs to fct and add to f_args
        cols = ()
        for key in col_keys:
            cols += (self[key],)
        f_args = cols + f_args

        # Define column and add it
        col = MaskedColumn(*args, data=fct(*f_args, **f_kwargs), **kwargs)
        self.add_column(col)
        
    def check_col_units(self, colname):
        
        col_units = self[colname].unit
        try:  # Check if col_units valid
            1. * Unit(col_units)
        except TypeError:
            print('Column has no units (unit = None)')
        except:  # Enter valid unit and refresh all table
            print("Column units '{}' are not".format(col_units) +
                          ' recognized by astropy.\n')
            print("Error message from astropy:")
            print_unit_error(str(col_units))
            print("-------------------------")
            gunit = input('***** Please enter the corresponding unit'
                          + ' recognized by astropy unit: ')
            self.correct_units(badunits=[str(col_units)], gunits=[gunit])
        
    def insert_value(self, pl_name, colname, value, reference, err=None, units=None, name_key=None,
                     extension='_refname', err_ext='err', ref_format='O', verbose=True):    
        """
        Insert new ``value`` in table for objet ``pl_name`` in
        column ``colname``. A ``reference`` is always needed.
        Errors (``err``) and ``units`` can be specified.

        Parameters
        ----------
        pl_name : str
            Planet name. The column for pl_name is specified by
            ``name_key`` or the table attribute ``main_col``.
            For NasaExoplanetArchive, it is ``pl_name``.
        colname : str
            Column name where to change the value
        value : same type as the dtype of the given Column
            Value to be inserted
        reference : str
            Reference name for the value. (Ex: 'Darveau-Bernier et al.')
        err : `None`, scalar or list of 2 values (plus/minus)
            Error for ``value``. If list, [err_plus, err_minus].
            If scalar, same error for both. Default is None. The name of 
            the error columns is given by ``colname``+``err_ext``+(1 or 2).
            If it doesn't match the actual err column name, it raises a
            warning and you should consider the method insert_error instead.
        units : str or astropy.units
            Units of the input ``value`` and ``err`` by extension. Optional.
            Default is the units of the given column.
        name_key : str or `None`
            Which column to be used to find ``pl_name``. Could be 'hd_name'
            for example if preferred. If None, default value is taken from
            Table attribute ``main_col``. An error will raise if none are defined.
        extension : str
            Extension for reference column. Default is '_refname'
        err_ext : str
            Extension for error columns. Default is 'err'.
        ref_format : numpy.dtype compatible value
            Will be given as an input to create reference column if it 
            does not exist yet. Correspond to dtype input in Column().
            Default is 'S15', a string format.

        Examples
        --------
        Example with bad error column name, units converted from rad to deg,
        specified name_key::
            >>> from transit_prediction.masterfile import MasterFile
            >>> data = MasterFile.read()
            >>> data.insert_value('HD 189733', 'ra', 4, 'Very-bad-reference et al.',
                                  err=0.1, units='rad', name_key='hd_name')
        Ouput::
            Setting column ra_refname to test
            Planet = b'HD 189733 b', column = ra
            Setting value to 229.1831180523293 deg
            Setting column [raerr1,raerr2] to [57.29577951308232 deg,57.29577951308232 deg]
            Last operation aborted, see WARNING
            WARNING: Wrong error column name. Did not write error value(s).
            Please use insert_error method and specify err colname. [warnings]

        """


        index = self.get_index(pl_name, name_key=name_key)[0]
        pl_name = self['pl_name'][index]
        log = []

        # Enter reference

        ref_col = colname + extension
        try:
            self[ref_col]
        except KeyError:
            self.add_column(MaskedColumn(length=len(self),
                                         mask=True,
                                         name=ref_col,
                                         dtype=ref_format))
        finally:
            log.append(
                'Setting column {} to {}'.format(ref_col, reference))
            if verbose:
                print(log[-1])
            self[ref_col][index] = reference

        # Check units and convert if needed
        
        self.check_col_units(colname)
        col_units = self[colname].unit

        if not units and col_units:
            if verbose:
                print('Default units are {}'.format(col_units))
                test = input('Does your input has the same units? (y/n):')
                if test != 'y':
                    raise UnitTypeError('Please enter the value in the good'+
                                       ' units or specify your units as input (units kwarg)')
            value = value * Unit(col_units)
            if err: err = err * Unit(col_units)
        elif units and col_units:
            value = (value * Unit(units)).to(col_units)
            if err: err = (err * Unit(units)).to(col_units)
        elif units and not col_units:
            raise UnitTypeError('Column has no units...' 
                                + 'see column description and leave units=None')

        # Enter new value
        
        log.append(
            'Planet = {}, column = {} \n'.format(pl_name, colname)
            + 'Setting value to {}'.format(value)
        )
        print(log[-1])
        try:
            self[colname][index] = value.value
        except AttributeError:
            self[colname][index] = value
            
        # Save log
        self.log.extend(log)
        log = []

        # Enter err value
        if err is not None:
            err_ext = colname + err_ext
            err_m_col, err_p_col = err_ext+'1', err_ext+'2'
            try:
                err_m, err_p = err
            except TypeError:
                err_m, err_p = err, err

            log.append(
                'Setting column [{},{}] to [{},{}]'\
                .format(err_m_col, err_p_col,err_m, err_p)
            )
            if verbose:
                print(log[-1])
            try:
                try:
                    self[err_m_col][index] = err_m.value
                    self[err_p_col][index] = err_p.value
                except AttributeError:
                    self[err_m_col][index] = err_m
                    self[err_p_col][index] = err_p
                # Save log
                self.log.extend(log)
            except KeyError:
                print('Last operation aborted, see WARNING')
                warnings.warn('Wrong error column name. Did not write error value(s).\n'+
                              'Please use insert_error method and specify err colname.',
                              AstropyUserWarning)

    def insert_error(self, pl_name, colname, err, units=None, name_key=None,
                     verbose=True):
        """
        Insert new ``err`` in table for objet ``pl_name`` in
        column ``colname``. A reference is not asked.
        ``units`` can be specified.

        Parameters
        ----------
        pl_name : str
            Planet name. The column for pl_name is specified by
            ``name_key`` or the table attribute ``main_col``.
            For NasaExoplanetArchive, it is ``pl_name``.
        colname : str
            Column name where to change the value
        err : same type as the dtype of the given Column
            Error value to be inserted
        units : str or astropy.units
            Units of the input ``err``. Optional.
            Default is the units of the given column.
        name_key : str or `None`
            Which column to be used to find ``pl_name``. Could be 'hd_name'
            for example if preferred. If `None`, default value is taken from
            Table attribute ``main_col``. An error will raise if none are defined.

        Examples
        --------
        None for now. Same usage as insert_value method but less complicated.

        """
        index = self.get_index(pl_name, name_key=name_key)[0]
        pl_name = self['pl_name'][index]

        # Check units and convert if needed
        
        self.check_col_units(colname)
        col_units = self[colname].unit
        if not units and col_units:
            if verbose:
                print('Default units are {}'.format(col_units))
                test = input('Does your input has the same units? (y/n):')
                if test != 'y':
                    raise UnitTypeError('Please enter the value in the good'+
                                       ' units or specify your units as input (units kwarg)')
            err = err * Unit(col_units)
        elif units and col_units:
            err = (err * Unit(units)).to(col_units)
        elif units and not col_units:
            raise UnitTypeError('Column has no units... see column description or leave units=None')

        # Enter new value

        self.log.append(
            'Planet = {}, column = {}\n'.format(pl_name, colname)
            + 'Setting value to {}'.format(err)
        )
        print(self.log[-1])
        try:
            self[colname][index] = err.value
        except AttributeError:
            self[colname][index] = err
            
    def add_column_beside(self, col, name_where, *args, where=1, **kwargs):
        
        index = self.index_column(name_where)
        index += where
        
        self.add_column(col, *args, index=index, **kwargs)


def difference(left, right):
    if isinstance(left, list) and isinstance(right, list):
        return list(set(left) - set(right))
    else:
        return NotImplemented


def intersection(self, other):
    if isinstance(self, list):
        return list(set(self).intersection(other))
    else:
        return NotImplemented
    
def print_unit_error(str_unit):
    
    try:
        Unit(str_unit)
    except ValueError as e:
        print(e)

                    
                 