# Python modules


# 3rd party modules


# Our modules



from vespa.common.util.import_ import Importer

import vespa.datasim.mrs_datasim as mrs_datasim


class DatasimImporter(Importer):
    def __init__(self, source):
        Importer.__init__(self, source, None, False)

    def go(self, add_history_comment=False):
        for element in self.root.getiterator("datasim"):
            self.found_count += 1
            datasim = mrs_datasim.Datasim(element)
            self.imported.append(datasim)

        self.post_import()
        
        return self.imported


class AnalysisHlsvdImporter(Importer):
    def __init__(self, source):
        Importer.__init__(self, source, None, False)

    def go(self, add_history_comment=False):
        for element in self.root.getiterator("analysis_export_hlsvd"):
            self.found_count += 1
            self.imported.append(element)

        self.post_import()

        return self.imported