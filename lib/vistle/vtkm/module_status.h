#ifndef VISTLE_VTKM_MODULE_STATUS_H
#define VISTLE_VTKM_MODULE_STATUS_H

#include <string>
#include <memory>

#include <vistle/core/messages.h>

#include "export.h"


//TODO: should this be moved to lib/vistle/module?

/**
 * @class ModuleStatus
 * @brief Base class for representing the status of a Vistle module.
 *
 * This class provides a common interface for representing the success or failure
 * of a module execution in classes which do not inherit from the Module class.
 * ModuleStatus objects can be passed to Vistle modules which then handle the status
 * further, e.g., by sending messages to the GUI or stopping their execution.
 */
class V_VTKM_EXPORT ModuleStatus {
protected:
    std::string msg;

public:
    ModuleStatus(const std::string &message);
    virtual ~ModuleStatus();

    // @return True if the module can continue its execution, false otherwise.
    virtual bool continueExecution() const = 0;

    // @return the message associated with the status
    const char *message() const;

    // @return The type of message to be sent to the GUI (e.g., Info, Warning, Error).
    virtual const vistle::message::SendText::TextType messageType() const = 0;
};

typedef std::unique_ptr<ModuleStatus> ModuleStatusPtr;

class V_VTKM_EXPORT SuccessStatus: public ModuleStatus {
public:
    SuccessStatus();

    bool continueExecution() const override;
    const vistle::message::SendText::TextType messageType() const override;
};

class V_VTKM_EXPORT ErrorStatus: public ModuleStatus {
public:
    ErrorStatus(const std::string &message);

    bool continueExecution() const override;
    const vistle::message::SendText::TextType messageType() const override;
};

class V_VTKM_EXPORT WarningStatus: public ModuleStatus {
public:
    WarningStatus(const std::string &message);

    bool continueExecution() const override;
    const vistle::message::SendText::TextType messageType() const override;
};

class V_VTKM_EXPORT InfoStatus: public ModuleStatus {
public:
    InfoStatus(const std::string &message);

    bool continueExecution() const override;
    const vistle::message::SendText::TextType messageType() const override;
};

/*
    @return A unique pointer to a ModuleStatus object representing a successful execution.
*/
ModuleStatusPtr V_VTKM_EXPORT Success();

/*
    @param message A C-string containing the error message.
    @return A unique pointer to a ModuleStatus object representing a failed execution.

*/
ModuleStatusPtr V_VTKM_EXPORT Error(const std::string &message);

/*
    @param message A C-string containing the warning message.
    @return A unique pointer to a ModuleStatus object representing a module that can continue its execution,
    but prints a warning message.
*/
ModuleStatusPtr V_VTKM_EXPORT Warning(const std::string &message);

/*
    @param message A C-string containing the informational message.
    @return A unique pointer to a ModuleStatus object representing a successful module that prints some 
    informational message.
*/
ModuleStatusPtr V_VTKM_EXPORT Info(const std::string &message);

#endif
