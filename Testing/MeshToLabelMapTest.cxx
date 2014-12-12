#include "itkTestMain.h"

#ifdef WIN32
#ifdef BUILD_SHARED_LIBS // variable declared in Cmake using "set_target_properties(... COMPILE_FLAGS -D...)"
#define MODULE_IMPORT __declspec(dllimport) // "dllimport" -> Dll = Dynamic library on windows -> DTIAtlasBuilderLib is static, so just set MODULE_IMPORT to nothing: Fails to link (linkage error) if dllimport because library is static and tries to link as dynamic
#else // BUILD_SHARED_LIBS
#define MODULE_IMPORT// define to empty
#endif // BUILD_SHARED_LIBS
#else // WIN32
#define MODULE_IMPORT// define to empty
#endif // WIN32

extern "C" MODULE_IMPORT int ModuleEntryPoint( int , char * [] ) ;
void RegisterTests()
{
  StringToTestFunctionMap["ModuleEntryPoint"] = ModuleEntryPoint ;
}
