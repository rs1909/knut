#include <QDialog>
#include <QTextEdit>
#include <QScrollBar>
#include <QVBoxLayout>
#include <QSettings>
#include <string>
#include <iostream>

class screenDialog : public QDialog
{
    Q_OBJECT
  public:
    screenDialog(QWidget *parent = 0) : QDialog(parent)
    {
      display = new QTextEdit;
      QVBoxLayout *layout = new QVBoxLayout;
      layout->addWidget(display);
      layout->setMargin(0);
      this->setLayout(layout);
      display->setLineWrapMode(QTextEdit::NoWrap);
      display->setReadOnly(true);
      QFont font("Courier");
      font.setFixedPitch(true);
      display->setFont(font);

      QSettings settings("Knut", "text window");
      QPoint pos = settings.value("pos", QPoint(200, 200)).toPoint();
      QSize size = settings.value("size", QSize(400, 400)).toSize();
      resize(size);
      move(pos);
    }
	~screenDialog()
	{
	  QSettings settings("Knut", "text window");
	  settings.setValue("pos", pos());
      settings.setValue("size", size());
	}
  public slots:
    void setTitle(const std::string& str)
    {
      display->setDocumentTitle(QString(str.c_str()));
    }
    void append(const std::string& str)
    {
      QTextCursor cur = display->textCursor ();
      cur.movePosition(QTextCursor::End);
      QScrollBar* bar = display->verticalScrollBar();
      const bool atEnd = (bar->maximum() == bar->value());
      cur.insertText(QString(str.c_str()));
      if (atEnd) bar->setValue(bar->maximum());
    }
    void setText(const QString& str)
    {
      display->setPlainText(str);
      display->setFontFamily("Courier");
    }
  private:
    QTextEdit *display;
};
