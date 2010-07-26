#include <QTextEdit>
#include <QScrollBar>
#include <QVBoxLayout>
#include <QSettings>
#include <QCloseEvent>
#include <string>
#include <iostream>

class screenDialog : public QTextEdit
{
    Q_OBJECT
  public:
    screenDialog(QWidget *parent = 0) : QTextEdit(parent)
    {
      setLineWrapMode(QTextEdit::NoWrap);
      setReadOnly(true);
      QFont font("Courier");
      font.setFixedPitch(true);
      setFont(font);

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
  signals:
    void windowClosed();
  protected:
  	void closeEvent(QCloseEvent *event)
 	{
 	  emit windowClosed();
      event->accept();
    }
  public slots:
    void setTitle(const std::string& str)
    {
      setDocumentTitle(QString(str.c_str()));
    }
    void append(const std::string& str)
    {
      QTextCursor cur = textCursor ();
      cur.movePosition(QTextCursor::End);
      QScrollBar* bar = verticalScrollBar();
      const bool atEnd = (bar->maximum() == bar->value());
      cur.insertText(QString(str.c_str()));
      if (atEnd) bar->setValue(bar->maximum());
    }
    void clearLastLine()
    {
      QTextCursor cur = textCursor ();
      cur.movePosition(QTextCursor::End);
      cur.select(QTextCursor::LineUnderCursor);
      cur.removeSelectedText();
    }
    void setText(const QString& str)
    {
      setPlainText(str);
      setFontFamily("Courier");
    }
};
