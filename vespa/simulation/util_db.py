# Python modules


# 3rd party modules

# Our modules
import vespa.common.util.db as common_util_db
from vespa.common.util.db import _Fetch


class Database(common_util_db.Database):
    """A database class that looks and works just like 
    common_util_db.Database. This exists to provide methods that are 
    specific to Simulation.
    
    In addition, this caches metabolites by default. (If you don't know
    what that means, don't worry about it.)
    """
    def __init__(self, filename, cache_metabolites=True):
        common_util_db.Database.__init__(self, filename, cache_metabolites)
        

    def append_simulations(self, experiment_id, simulations):
        """Given an experiment id and an iterable of mrs_simulation objects,
        attaches those simulations to the experiment.
        
        This function does not use a transaction. The caller should.
        """
        self.begin_transaction()
        self._insert_simulations(experiment_id, simulations)
        self.commit()
        

    def delete_experiments(self, experiment_ids, use_transaction=True):
        """Given a single experiment id or a list of them, deletes the 
        experiment(s) and associated simulations.

        If use_transaction=True, the database will wrap each set of deletes 
        in a transaction. There's no reason to set this to False unless you're
        calling this from code which already has a transaction open. 
        """
        if isinstance(experiment_ids, str):
            experiment_ids = [experiment_ids]
                
        for id_ in experiment_ids:
            if use_transaction:
                self.begin_transaction()
            
            # Clean up in reverse order -- first simulations, then dims, 
            # then params, then experiment_metabs, and lastly the experiment 
            # itself. This preserves referential integrity as I work which 
            # isn't necessary but I like to do it anyway.

            # Delete simulations
            sql = """DELETE FROM
                        simulations
                     WHERE
                        simulations.dims_id IN
                            (SELECT
                                id
                             FROM
                                experiment_dims
                             WHERE
                                experiment_id = ?)
                   """
            self._execute(sql, id_, _Fetch.NONE)
            
            # Delete dims
            sql = """DELETE FROM
                        experiment_dims
                     WHERE
                        experiment_id = ?
                   """
            self._execute(sql, id_, _Fetch.NONE)
            
            # Delete parameters
            sql = """DELETE FROM
                        experiment_user_static_parameters
                     WHERE
                        experiment_id = ?
                   """
            self._execute(sql, id_, _Fetch.NONE)
            
            # Delete metab references
            sql = """DELETE FROM
                        experiment_metabolites
                     WHERE
                        experiment_id = ?
                   """
            self._execute(sql, id_, _Fetch.NONE)
            
            # Finally, delete the experiment
            sql = """DELETE FROM
                        experiments
                     WHERE
                        id = ?
                   """
            self._execute(sql, id_, _Fetch.NONE)

            if use_transaction:
                self.commit()
        

    def delete_metabolites(self, metabolite_ids):
        """Given a single metabolite id or a list of them, deletes the 
        metabolite(s) and associated spins & J couplings. 
        
        If any of the metabolites are in use by experiments, the deletion
        will fail gracelessly.
        """
        if isinstance(metabolite_ids, str):
            metabolite_ids = [metabolite_ids]
                
        for id_ in metabolite_ids:
            self.begin_transaction()
            
            # Failsafe here -- metabolites that are referenced by an 
            # experiment must not be deleted and should never be passed to
            # this function.
            sql = """SELECT
                        count(*)
                     FROM
                        experiment_metabolites
                     WHERE
                        metabolite_id = ?"""
            usage_count = self._execute(sql, id_, _Fetch.SINGLETON)
            
            assert(not usage_count)

            # Clean up in reverse order -- first J Couplings, then the
            # spins to which they refer, and lastly the metab itself. This
            # preserves referential integrity as I work which isn't necessary
            # but I like to do it anyway.

            # Find the spins referenced by this metabolite
            sql = """SELECT
                        id
                     FROM
                        metabolite_spins
                     WHERE
                        metabolite_id = ?"""
            spin_ids = self._execute(sql, id_, _Fetch.COLUMN)

            # Delete j couplings
            placeholders = common_util_db._build_placeholders(len(spin_ids))
            sql = """DELETE FROM
                        j_couplings
                     WHERE
                        spin1_id IN (%s) OR
                        spin2_id IN (%s)
                   """ % (placeholders, placeholders)

            self._execute(sql, spin_ids + spin_ids, _Fetch.NONE)

            # Delete spins
            sql = """DELETE FROM
                        metabolite_spins
                     WHERE
                        id IN (%s)
                   """ % placeholders
            self._execute(sql, spin_ids, _Fetch.NONE)


            # Delete metab
            sql = """DELETE FROM
                        metabolites
                     WHERE
                        id = ?
                   """
            self._execute(sql, id_, _Fetch.NONE)

            self.commit()
        
            if id_ in self._metabolite_cache:
                del self._metabolite_cache[id_]


    def delete_pulse_sequences(self, pulse_sequence_ids):
        """Given a single pulse sequence id or a list of them, deletes the 
        pulse_sequence(s) and associated params. 
        
        If any of the pulse sequences are in use by experiments, the deletion
        will fail. 
        """
        if isinstance(pulse_sequence_ids, str):
            pulse_sequence_ids = [pulse_sequence_ids]
                
        for id_ in pulse_sequence_ids:
            self.begin_transaction()
            
            # Failsafe here -- pulse sequences that are referenced by an 
            # experiment must not be deleted and should never be passed to
            # this function.
            sql = """SELECT
                        count(*)
                     FROM
                        experiments
                     WHERE
                        pulse_sequence_id = ?"""
            usage_count = self._execute(sql, id_, _Fetch.SINGLETON)
            
            assert(not usage_count)

            # Delete loops, then params, then the pulse sequence.
            sql = """DELETE FROM
                        pulse_sequence_loops
                     WHERE
                        pulse_sequence_id = ?"""
            self._execute(sql, id_, _Fetch.NONE)
        
            sql = """DELETE FROM
                        pulse_sequence_user_static_parameters
                     WHERE
                        pulse_sequence_id = ?"""
            self._execute(sql, id_, _Fetch.NONE)

            sql = """DELETE FROM
                        pulse_sequence_pulse_designs
                     WHERE
                        pulse_sequence_id = ?"""
            self._execute(sql, id_, _Fetch.NONE)

            sql = """DELETE FROM
                        pulse_sequences
                     WHERE
                        id = ?
                   """
            self._execute(sql, id_, _Fetch.NONE)

            self.commit()


    def replace_experiment(self, experiment):
        """
        Deletes and re-inserts an experiment with the same UUID. The 
        experiment must exist in the database.

        This is as close as we get to update_experiment().
        """
        self.begin_transaction()

        # Ensure the experiment exists.
        assert(bool(self.count_experiments(experiment.id)))
        
        self.delete_experiments(experiment.id, False)

        self.insert_experiment(experiment, False)
        
        self.commit()


    def update_metabolite(self, metabolite):
        """Saves the metabolite described by the metabolite param. The 
        metabolite must already exist in the database.
        """
        self.begin_transaction()
        
        # There's two ways I could do an update of a metab's spins & 
        # J couplings. One would be to compare each spin & J coupling on the
        # metabolite param object with what's in the database and resolve
        # the differences. The second method is a lot easier -- assume
        # they've changed, delete them from the database and write them anew.
        
        # Delete the spins & J couplings currently associated with this 
        # metab in the database.
        sql = """SELECT
                    id
                 FROM
                    metabolite_spins
                 WHERE
                    metabolite_id = ?"""
        spin_ids = self._execute(sql, metabolite.id, _Fetch.COLUMN)
        
        if spin_ids:
            placeholders = common_util_db._build_placeholders(len(spin_ids))
            sql = """DELETE FROM
                        j_couplings
                     WHERE
                        spin1_id in (%s) OR
                        spin2_id in (%s)
                  """ % (placeholders, placeholders)
            self._execute(sql, spin_ids + spin_ids, _Fetch.NONE)

            sql = """DELETE FROM
                        metabolite_spins
                     WHERE
                        id in (%s)
                  """ % placeholders
            self._execute(sql, spin_ids, _Fetch.NONE)
                  
        # Add spins
        sql = """INSERT INTO
                    metabolite_spins (metabolite_id, isotope, chemical_shift,
                                      display_order)
                 VALUES
                    (?, ?, ?, ?)
              """
        for i, spin in enumerate(metabolite.spins):
            params = (metabolite.id, spin.isotope, spin.chemical_shift, i)

            spin.id = self._execute(sql, params, _Fetch.LAST_ID)

        # Add J couplings
        sql = """INSERT INTO
                    j_couplings (value, spin1_id, spin2_id)
                 VALUES
                    (?, ?, ?)
              """
        for j_coupling in metabolite.j_couplings:
            params = (j_coupling.value, j_coupling.spin1.id, j_coupling.spin2.id)

            j_coupling.id = self._execute(sql, params, _Fetch.LAST_ID)
        
        # Update the metabolite
        sql = """UPDATE
                    metabolites
                 SET
                    name = ?, is_public = ?, creator = ?, comment = ?, 
                    deactivated = ?
                 WHERE
                    id = ?
              """
        self._execute(sql, (metabolite.name, metabolite.is_public,
                            metabolite.creator, metabolite.comment,
                            metabolite.deactivated, metabolite.id),
                      _Fetch.NONE)
                      
        self.commit()

        if self.cache_metabolites:
            # Ensure that the cache is up-to-date
            self._update_metabolite_cache(metabolite.id)

        
    def update_pulse_sequence(self, pulse_sequence):
        """Saves the pulse_sequence described by the param. The 
        pulse sequence must already exist in the database.
        """
        self.begin_transaction()
        
        # Delete the existing loops
        sql = """DELETE FROM
                    pulse_sequence_loops
                 WHERE
                    pulse_sequence_id = ?"""
        self._execute(sql, pulse_sequence.id, _Fetch.NONE)
        
        # Insert the loops      
        sql = """INSERT INTO
                    pulse_sequence_loops
                        (pulse_sequence_id, label, display_order)
                 VALUES
                    (?, ?, ?)
              """
        for i, label in enumerate(pulse_sequence.loop_labels):
            if label:
                params = (pulse_sequence.id, label, i) 

            self._execute(sql, params, _Fetch.NONE)

        # Delete the existing params
        sql = """DELETE FROM
                    pulse_sequence_user_static_parameters
                 WHERE
                    pulse_sequence_id = ?"""
        self._execute(sql, pulse_sequence.id, _Fetch.NONE)
        
        # Add new ones
        sql = """INSERT INTO
                    pulse_sequence_user_static_parameters 
                        (pulse_sequence_id, type, name, default_value,
                         display_order)
                 VALUES
                    (?, ?, ?, ?, ?)
              """
        for i, parameter in enumerate(pulse_sequence.user_static_parameters):
            params = (pulse_sequence.id, parameter.type, parameter.name, 
                      parameter.default, i)

            parameter.id = self._execute(sql, params, _Fetch.NONE)

        # Delete the existing pulse design references
        sql = """DELETE FROM
                    pulse_sequence_pulse_designs
                 WHERE
                    pulse_sequence_id = ?"""
        self._execute(sql, pulse_sequence.id, _Fetch.NONE)
        
        # Add new pulse project references
        sql = """INSERT INTO
                    pulse_sequence_pulse_designs
                        (pulse_sequence_id, pulse_design_id, progression)
                 VALUES
                    (?, ?, ?)
              """
        for i, pulse_design in enumerate(pulse_sequence.pulse_projects):
            self._execute(sql, (pulse_sequence.id, pulse_design.id, i))



        # Update the pulse_sequence itself
        sql = """UPDATE
                    pulse_sequences
                 SET
                    name = ?, comment = ?, is_public = ?,
                    sequence_code = ?, binning_code = ?
                 WHERE
                    id = ?
              """
        params = (pulse_sequence.name, pulse_sequence.comment,
                  pulse_sequence.is_public,
                  pulse_sequence.sequence_code, pulse_sequence.binning_code, 
                  pulse_sequence.id) 
                 
        self._execute(sql, params, _Fetch.NONE)
                      
        self.commit()
        
