# Python modules


# 3rd party modules


# Our modules
import vespa.common.util.db as common_db
import vespa.common.util.misc as util_misc
from vespa.common.util.db import _Fetch
import vespa.common.rfp_transform_kernel as rfp_transform_kernel
import vespa.common.rfp_machine_specs as rfp_machine_specs
import vespa.common.constants as constants




class Database(common_db.Database):
    """A database class that looks and works just like 
    vespa.common.util.db.Database. This exists to provide methods that are 
    specific to Pulse.
    """
    def __init__(self, filename):
        common_db.Database.__init__(self, filename)


##############                       ###########################  
##############  Pulse Design Methods ###########################
##############                       ###########################

    def exists_pulse_design_in_database(self, pdid):
        """ Lets first make sure the pulse_design is there! """
        sql = """SELECT
                    *
                 FROM
                    pulse_designs
                 WHERE
                    id = ?
              """
        return bool(self._execute(sql, pdid, fetch=_Fetch.ONE))    


    def delete_pulse_design(self, pulse_design_id, use_transactions=True):
        """Delete the pulse_design with the given pulse_design_id."""
        
        if use_transactions:
            self.begin_transaction()
            
        sql = """SELECT
                    *
                 FROM
                    pulse_designs
                 WHERE
                    id = ?
              """
        row = self._execute(sql, pulse_design_id, fetch=_Fetch.ONE)
  
        self.delete_machine_specs(row["machine_specs_id"]) 

        sql = """SELECT
                    *
                 FROM
                    transforms
                 WHERE
                    pulse_design_id = ?
              """
        rows = self._execute(sql, (pulse_design_id), fetch=_Fetch.ALL)

        for row in rows:
            self._delete_transform_parameters(row['id'])
              
        for row in rows:
            self._delete_rf_result(row["rf_result_id"])
              
        # Delete Transforms
        sql = """DELETE FROM
                    transforms
                 WHERE
                    pulse_design_id = ?
              """
        self._execute(sql, (pulse_design_id), fetch=_Fetch.NONE)

        sql = """DELETE FROM
                    pulse_designs
                 WHERE
                    id = ?
               """                       
        self._execute(sql, pulse_design_id, _Fetch.NONE)

        if use_transactions:
            self.commit()
    
    
    def delete_pulse_designs(self, pulse_design_ids):
        self.begin_transaction()
        for pdid in pulse_design_ids:
            self.delete_pulse_design(pdid, False)
        self.commit()


    def replace_pulse_design(self, pulse_design, results_to_skip=[ ]):
        """
        Deletes and re-inserts a pulse_design with the same UUID. 
        The pulse design must exist in the database.

        This is as close as we get to update_pulse_design().

        See insert_pulse_design() for an explanation of results_to_skip
        """
        
        self.begin_transaction()

        assert(self.exists_pulse_design_in_database(pulse_design.id))
        
        self.delete_pulse_design(pulse_design.id, False)

        self.insert_pulse_design(pulse_design, False, results_to_skip)
        
        self.commit()






##############                               ###########################  
##############  Pulse Design Preview Methods ###########################
##############                               ###########################
    
# all are in the main db.py because they get used in /common modules

        
  
##############                    ###########################  
##############  Transform Methods ###########################
##############                    ###########################

#   Note.  delete_transforms() and fetch_transforms() are not needed as
#          they are integral to the xxx_pulse_design methods

    def _delete_transform_parameters(self, transform_id):
        """
        Currently this code is only called from within a transaction, so
        it doesn't need to open a transaction itself.

        Given a single transform id, delete the transform parameters 
        associated with it. 
        
        """
        sql = """DELETE FROM
                    transform_parameters
                 WHERE
                    transform_id = ?"""
        self._execute(sql, transform_id, _Fetch.NONE)






##############                           ###########################  
##############  Transform Kernel Methods ###########################
##############                           ###########################

    def get_transform_kernel_menu_items(self):
        """ get labels, and ids for Transform menu """
        sql = """SELECT
                    id, menu_label
                 FROM
                    transform_kernels
                 WHERE
                    1 = 1
                 ORDER BY 
                    menu_label
              """
        rows = self._execute(sql, None)
        
        rows = [(item['id'], item['menu_label']) for item in rows]
        
        return rows


    def get_transform_kernel_create_menu_items(self):
        """ get labels, and ids for Transform menu """
        sql = """SELECT
                    id, menu_label
                 FROM
                    transform_kernels
                 WHERE
                    type = 'Create Transform'
                 ORDER BY 
                    menu_label
              """
        rows = self._execute(sql, None)
        
        rows = [(item['id'], item['menu_label']) for item in rows]
        
        return rows


    def get_transform_kernel_modify_menu_items(self):
        """ get labels, and ids for Transform menu """
        sql = """SELECT
                    id, type, menu_label
                 FROM
                    transform_kernels
                 WHERE
                    type != 'Create Transform'
                 ORDER BY 
                    menu_label
              """
        rows = self._execute(sql, None)
        
        rows = [(item['id'], item['menu_label']) for item in rows]
        
        return rows
    

    def exists_transform_kernel_in_database(self, tkid):
        """ Lets first make sure the transform_kernel is there! """
        sql = """SELECT
                    *
                 FROM
                    transform_kernels
                 WHERE
                    id = ?
              """
        return bool(self._execute(sql, tkid, fetch=_Fetch.ONE))    


    def delete_transform_kernels(self, transform_kernel_ids):
        """
        Given a single transform_kernel id or a list of them, deletes the 
        transform_kernel(s) and associated params. 
        
        If any of the pulse sequences are in use by experiments, the deletion
        will fail. 
        """
        if isinstance(transform_kernel_ids, str):
            transform_kernel_ids = [transform_kernel_ids]
                
        for id_ in transform_kernel_ids:
            self.begin_transaction()
            
            # Failsafe here -- transform_kernel_ids that are referenced by a 
            # transform must not be deleted and should never be passed to
            # this function.
            sql = """SELECT
                        count(*)
                     FROM
                        transforms
                     WHERE
                        transform_kernel_id = ?"""
            usage_count = self._execute(sql, id_, _Fetch.SINGLETON)
             
            #assert(not usage_count)

            # Delete transform_kernel_controls, then the transform_kernel.
            sql = """DELETE FROM
                        transform_kernel_controls
                     WHERE
                        transform_kernel_id = ?"""
            self._execute(sql, id_, _Fetch.NONE)

            sql = """DELETE FROM
                        transform_kernels
                     WHERE
                        id = ?
                   """
            self._execute(sql, id_, _Fetch.NONE)

            self.commit()


    def replace_transform_kernel(self, transform_kernel):
        """
        Deletes and re-inserts a transform_kernel with the same UUID. 
        The transform_kernel must exist in the database.

        This is as close as we get to update_transform_kernel().

        """
        self.begin_transaction()

        assert(self.exists_transform_kernel_in_database(transform_kernel.id))
        
        self.delete_transform_kernels(transform_kernel.id)

        self.insert_transform_kernel(transform_kernel)
        
        self.commit()
                 
    
    
    def fetch_transform_kernels(self):
        """
        Given a transform kernel id, returns the associated 
        transform_kernel object.
        
        """
        sql = """SELECT
                    transform_kernels.*
                 FROM
                    transform_kernels
                 WHERE
                    1 = 1
                 ORDER BY 
                    name
              """
        rows = self._execute(sql, None)

        # Will inflate row into a transform_kernel.
        transform_kernels = [rfp_transform_kernel.TransformKernel(row) for row in rows]
        
        for transform_kernel in transform_kernels:
            # Fetch pulse sequence params. (There might not be any.)
            sql = """SELECT
                        *
                     FROM
                        transform_kernel_controls
                     WHERE
                        transform_kernel_id = ?
                     ORDER BY
                        display_order
                     """
            rows = self._execute(sql, transform_kernel.id)

            # Some special attention is required here. Unlike most of our
            # tables/objects, there's not a 1:1 map between the column and 
            # attribute names for pulse sequence params. The table 
            # transform_kernel_controls contains a column called 
            # default_ while the TransformKernelControl class has an 
            # attribute called default. Similarly for name_, type_ and
            # variable_ columns. We can't use the latter as a column 
            # name because it's a reserved word in SQL, and neither can we 
            # use it as a column name alias. 
            # So here I rename the column in the returned row. Remember 
            # that the rows returned by self._execute() are _BetterRow objects
            # and need to be turned into proper Python dicts before I 
            # manipulate them.
            rows = [dict(row) for row in rows]
            
            for row in rows:
                row["name"]     = row["name_"]
                row["type"]     = row["type_"]
                row["default"]  = row["default_"]
                row["variable"] = row["variable_"]
                del row["name_"]
                del row["type_"]
                del row["default_"]
                del row["variable_"]

            transform_kernel.transform_kernel_controls = [rfp_transform_kernel.TransformKernelControl(row) for row in rows]

            transform_kernel.referrers = self.fetch_transform_kernel_referrers(transform_kernel.id)

        return transform_kernels



##############                     ###########################  
##############  Rf Results Methods ###########################
##############                     ###########################

    def _delete_rf_result(self, rf_result_id):    
        # Currently this code is only called from within a transaction, so
        # it doesn't need to open a transaction itself.
        if not rf_result_id:
            return
 
        sql = """SELECT
                    *
                 FROM
                    rf_results
                 WHERE
                    id = ?
              """
        row = self._execute(sql, rf_result_id, fetch=_Fetch.ONE) 
        
        self._delete_opcon_state(row["opcon_state_id"]) 
 
        # Now delete the "container" table for this result.
        sql = '''DELETE FROM
                    rf_results
                 WHERE
                    id = ?
              '''   
        self._execute(sql, rf_result_id, _Fetch.NONE)  


    
        
        
##############                        ###########################  
##############  Machine Specs Methods ###########################
##############                        ###########################

    def delete_machine_specs(self, ids):
        """
        Given a single machine spec id or a list of them, deletes the 
        machine specs. It doesn't matter if any/all of the machine specs
        are templates.
        """
        # Make ids a list if it isn't one already
        if not util_misc.is_iterable(ids, False):
            ids = (ids, )
            
        sql = """DELETE FROM
                    machine_specs
                 WHERE
                    id IN (%s)
               """ % common_db._build_placeholders(len(ids))

        self._execute(sql, ids, _Fetch.NONE)

        
        
    def fetch_machine_specs(self, id_):
        """
        Given the id of a machine specs object, returns that object. It
        doesn't matter if it's a template or not.
        """
        
        sql = """SELECT
                    *
                 FROM
                    machine_specs
                 WHERE
                    id = ?
              """                                      
        row = self._execute(sql, id_, _Fetch.ONE)
        
        if row["is_template"]:
            return rfp_machine_specs.MachineSpecsTemplate(row)
        else:
            return rfp_machine_specs.MachineSpecs(row)


    def fetch_machine_specs_templates(self):
        """Returns a list of machine_specs templates ordered by name"""
        sql = """SELECT
                    *
                 FROM
                    machine_specs
                 WHERE
                    is_template = 1
                 ORDER BY
                    name
              """                                      
        rows = self._execute(sql)
        
        return [rfp_machine_specs.MachineSpecsTemplate(row) for row in rows]



    def update_machine_specs(self, machine_specs):
        """
        Saves the machine specs described by the input parameter. 
        The machine specs must already exist in the database.
        
        It doesn't matter if the machine specs object is a template or not.
        
        """
        # Templates are a little different from the specs associated with
        # a pulse design.
        is_template = isinstance(machine_specs, rfp_machine_specs.MachineSpecsTemplate)

        name = machine_specs.name if is_template else None
        is_default = machine_specs.is_default if is_template else False
        # I translate machine type to the correct constant, but only if it's
        # not freeform text.
        machine_type = machine_specs.machine_type
        if machine_type in constants.MachineType.ALL:
            machine_type = machine_type['db']
        
        sql = """UPDATE
                    machine_specs
                 SET
                    name = ?,
                    is_default = ?,
                    machine_type = ?, 
                    field_strength = ?, 
                    max_b1_field = ?,
                    zero_padding = ?, 
                    min_dwell_time = ?,
                    dwell_time_increment = ?, 
                    gradient_raster_time = ?,
                    gradient_slew_rate = ?,
                    gradient_maximum = ?
                 WHERE
                    id = ?
              """
        
        sql_params = (name, 
                      is_default, 
                      machine_type,
                      machine_specs.field_strength,
                      machine_specs.max_b1_field,
                      machine_specs.zero_padding,
                      machine_specs.min_dwell_time,
                      machine_specs.dwell_time_increment,
                      machine_specs.gradient_raster_time,
                      machine_specs.gradient_slew_rate,
                      machine_specs.gradient_maximum,
                      machine_specs.id)

        self._execute(sql, sql_params, _Fetch.NONE)



##############                       ###########################  
##############  Opcon States Methods ###########################
##############                       ###########################

    def _delete_opcon_parameters(self, parameters_id):

        sql = """DELETE FROM
                    opcon_parameters
                 WHERE
                    id = ?
              """
        self._execute(sql, parameters_id, _Fetch.NONE)


    def _delete_opcon_state(self, opcon_state_id):
        
        if not opcon_state_id:
            return
 
        # Delete all the individual data points
        sql = '''DELETE FROM
                    deltab1_points
                 WHERE
                    opcon_state_id = ?
              '''          
        self._execute(sql, opcon_state_id, _Fetch.NONE)         
               
        # Delete the residual error history
        sql = '''DELETE FROM
                    opcon_residual_errors
                 WHERE
                    opcon_state_id = ?
              '''          
        self._execute(sql, opcon_state_id, _Fetch.NONE)

        # Delete the OC state itself
        sql = '''DELETE FROM
                    opcon_states
                 WHERE
                    id = ?                        
              '''        
        self._execute(sql, opcon_state_id, _Fetch.NONE)  


