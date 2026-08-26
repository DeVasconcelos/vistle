#include <cstring>
#include <exception>
#include <string>
#include <sstream>

#include <viskores/cont/Error.h>
#include <viskores/cont/ErrorBadAllocation.h>
#include <viskores/cont/ErrorBadDevice.h>
#include <viskores/cont/ErrorBadType.h>
#include <viskores/cont/ErrorBadValue.h>
#include <viskores/cont/ErrorExecution.h>
#include <viskores/cont/ErrorFilterExecution.h>
#include <viskores/cont/ErrorInternal.h>

#include "vtkm_module_utils.h"

namespace vistle {
namespace utils {

bool isValid(const vistle::Module &module, const ModuleStatusPtr &status)
{
    if (strcmp(status->message(), ""))
        module.sendText(status->messageType(), status->message());

    return status->continueExecution();
}

bool tryToExecuteFilter(const vistle::Module &module, const std::unique_ptr<viskores::filter::Filter> &filter,
                        const viskores::cont::DataSet &inputDataset, viskores::cont::DataSet &outputDataset)
{
    std::string kind, description, message, backtrace;

    try {
        try {
            outputDataset = filter->Execute(inputDataset);
            return true;
        } catch (const viskores::cont::ErrorBadAllocation &error) {
            kind = "memory allocation error";
            description = "A memory allocation error occurred while executing the filter";
            throw;
        } catch (const viskores::cont::ErrorBadDevice &error) {
            kind = "operation not supported by execution device";
            description = "The filter attempted to perform an operation that is not supported by the execution device";
            throw;
        } catch (const viskores::cont::ErrorBadType &error) {
            kind = "unsupported data type";
            description = "An unsupported data type was encountered while executing the filter";
            throw;
        } catch (const viskores::cont::ErrorBadValue &error) {
            kind = "unsupported data value";
            description = "An invalid value was encountered while executing the filter";
            throw;
        } catch (const viskores::cont::ErrorExecution &error) {
            kind = "execution environment error";
            description = "An error occurred in the execution environment while executing the filter";
            throw;
        } catch (const viskores::cont::ErrorFilterExecution &error) {
            kind = "filter setup error";
            description = "The filter has not been set up correctly";
            throw;
        } catch (const viskores::cont::ErrorInternal &error) {
            kind = "internal error";
            description = "An internal error occurred while executing the filter, indicating a bug in Viskores";
            throw;
        }
    } catch (const viskores::cont::Error &error) {
        if (kind.empty())
            kind = "unknown error";
        kind = "Viskores " + kind;
        message = error.GetMessage();
        backtrace = error.GetStackTrace();
    } catch (std::exception &error) {
        kind = "standard exception ";
        kind += typeid(error).name();
        message = error.what();
    } catch (...) {
        kind = "unknown exception";
    }

    std::stringstream msg;
    msg << "Execution of a Viskores filter failed with a " << kind << " exception";
    if (!description.empty())
        msg << " (" << description << ")";
    if (!message.empty())
        msg << ": " << message;
    module.sendError(msg.str());

    if (!backtrace.empty()) {
        msg << "\nBacktrace:\n" << backtrace;
    }

    std::cerr << msg.str() << std::endl;

    return false;
}

} // namespace utils
} // namespace vistle
