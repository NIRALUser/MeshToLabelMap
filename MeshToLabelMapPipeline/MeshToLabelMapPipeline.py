from slicer.ScriptedLoadableModule import ScriptedLoadableModule

try:
  from PipelineModulesLib.CLIModuleWrapping import PipelineCLI
  PipelineCLI("MeshToLabelMap", inputArgName="mesh", excludeArgs=['reference'])
except ImportError:
  pass

#
# MeshToLabelMapPipeline
#

class MeshToLabelMapPipeline(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "MeshToLabelMapPipeline"
    self.parent.categories = ["Surface Models/Pipelines"]
    self.parent.dependencies = ["MeshToLabelMap"]
    self.parent.contributors = ["Connor Bowley (Kitware, Inc)"]
    self.parent.helpText = "This module exists to provide a pipeline for mesh to labelmap"
    self.parent.acknowledgementText = ""
    self.parent.hidden = True
