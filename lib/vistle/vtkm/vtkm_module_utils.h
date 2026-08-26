#ifndef VISTLE_VTKM_VTKM_MODULE_UTILS_H
#define VISTLE_VTKM_VTKM_MODULE_UTILS_H

#include <memory> // for unique_ptr

#include <viskores/filter/Filter.h>

#include <vistle/module/module.h>

#include "module_status.h"

namespace vistle {
namespace utils {
/*
    @brief Checks if the module status is valid and sends a message to the GUI if necessary.

    This function checks if the provided ModuleStatus object indicates that the module can
    continue its execution. If the status contains a message, it sends the message to the 
    GUI using the Module's send-methods.

    @return True if the module can continue its execution, false otherwise.
*/
bool V_VTKM_EXPORT isValid(const vistle::Module &module, const ModuleStatusPtr &status);

/*
    @brief Attempts to execute a Viskores filter and handles any exceptions that may occur.

    This function attempts to execute the provided Viskores filter on the input dataset. If
    an exception occurs during execution, it captures the exception details and sends an 
    appropriate message to the GUI using the Module's send-methods.

    @return True if the filter was executed successfully, false otherwise.
*/
bool V_VTKM_EXPORT tryToExecuteFilter(const vistle::Module &module,
                                      const std::unique_ptr<viskores::filter::Filter> &filter,
                                      const viskores::cont::DataSet &inputDataset,
                                      viskores::cont::DataSet &outputDataset);
} // namespace utils
} // namespace vistle

#endif // VISTLE_VTKM_VTKM_MODULE_UTILS_H
