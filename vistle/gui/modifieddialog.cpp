#include "modifieddialog.h"

#include <QPushButton>
#include <QDialogButtonBox>

#include <cassert>

namespace gui {

ModifiedDialog::ModifiedDialog(const QString &reason, QWidget *parent)
: QMessageBox(parent)
{
   setWindowModality(Qt::WindowModal);

   setIcon(QMessageBox::Question);

   setWindowTitle(QString("%1").arg(reason));

   setText("The data flow network has been modified.");
   if (reason == "Open") {
      setInformativeText("You will lose your changes, if you load another data flow network without saving.");
   } else if (reason == "New") {
      setInformativeText("You will lose your changes, if you clear the data flow network without saving.");
   } else if (reason == "Quit") {
      setInformativeText("You will lose your changes, if you quit without saving.");
   } else {
      setInformativeText("You will lose your changes, if you continue without saving.");
   }

   setStandardButtons(QMessageBox::Cancel | QMessageBox::Save);
   m_discardButton = addButton(reason, QMessageBox::DestructiveRole);
   setDefaultButton(m_discardButton);
}

int ModifiedDialog::exec() {

   int result = QMessageBox::exec();
   if (result == QMessageBox::Cancel)
      return result;
   if (result == QMessageBox::Save)
      return result;

   if (clickedButton() == m_discardButton)
      return QMessageBox::Discard;

   assert("QMessageBox result not handled" == 0);
   return QMessageBox::Cancel;
}

} // namespace gui
