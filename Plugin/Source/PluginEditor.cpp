#include "PluginProcessor.h"
#include "PluginEditor.h"
LiveProgEditor::LiveProgEditor(LiveProgrammableDSP &p) : AudioProcessorEditor(&p), processor(p), table(processor.ptrVM_var_name, processor.ptrVM_var_value, processor.ptrNumBlocks)
{
	fileComp = new FilenameComponent("fileComp",
		{},                       // current file
		false,                    // can edit file name,
		false,                    // is directory,
		false,                    // is for saving,
		{},                       // browser wildcard suffix,
		{},                       // enforced suffix,
		"Select EEL/Lua routine to open");  // text when nothing selected
	addAndMakeVisible(fileComp);
	fileComp->addListener(this);
	fileComp->setDefaultBrowseTarget(File::getSpecialLocation(File::userDocumentsDirectory));

	textContent.reset(new Pop_Texteditor());
	addAndMakeVisible(textContent.get());
	textContent->setMultiLine(true, false);
	textContent->setReadOnly(false);
	textContent->setCaretVisible(true);
	textContent->setReturnKeyStartsNewLine(true);
	textContent->setFont(Font(17));
	textContent->fileCompPtr = fileComp;


	consoleOutput.reset(new Pop_Texteditor());
	addAndMakeVisible(consoleOutput.get());
	consoleOutput->setMultiLine(true, false);
	consoleOutput->setReadOnly(false);
	consoleOutput->setCaretVisible(true);
	consoleOutput->setReturnKeyStartsNewLine(true);
	consoleOutput->setFont(Font(17));
	consoleOutput->fileCompPtr = fileComp;

	addAndMakeVisible(&playButton);
	playButton.setButtonText("Compile & Run");
	playButton.onClick = [this]
	{
		table.allowUpdate = 0;
		processor.restoredText = textContent->getText();
		processor.playButtonClicked(processor.restoredText);
		table.ptrVar_name = processor.ptrVM_var_name;
		table.ptrVar_value = processor.ptrVM_var_value;
		table.ptrNB = processor.ptrNumBlocks;
//		errDisp->setText(processor.errorMsg, NotificationType::dontSendNotification);
		table.allowUpdate = 1;
		table.timerCallback();
		table.repaint();
	};
	playButton.setColour(TextButton::buttonColourId, Colours::green);
	playButton.setEnabled(true);
//	errDisp = new Label("errLbl", "");
//	errDisp->setEnabled(true);
//	addAndMakeVisible(errDisp);
	setSize(600, 400);
	jtooltipWindow->setMillisecondsBeforeTipAppears(1000);
	//addAndMakeVisible(jsocialButtons);
	if (processor.restoredText.isNotEmpty())
		textContent->setText(processor.restoredText);
//	if (processor.errorMsg.isNotEmpty())
//		errDisp->setText(processor.errorMsg, NotificationType::dontSendNotification);

	setResizable(true, true);
	setResizeLimits(400, 350, 800, 750);
	addAndMakeVisible(table);
	// set plugin's UI window size
	setSize(500, 450);
	startTimer(UPDATEINTERVAL);
}
LiveProgEditor::~LiveProgEditor()
{
	delete fileComp;
//	delete errDisp;
}
#ifdef CUSTOM_CMD
extern "C"
{
	int needUpdate = 0;
	int circularStrCurrent = 0;
	char history[HISTORY_COUNT][MAX_CMD_LEN] = { 0 };
	void writeCircularStringBuf(char *cmdCur)
	{
		memset(history[circularStrCurrent], 0, MAX_CMD_LEN);
		memcpy(history[circularStrCurrent], cmdCur, strlen(cmdCur));
		circularStrCurrent = (circularStrCurrent + 1) % HISTORY_COUNT;
		needUpdate = 1;
	}
}
#endif
void LiveProgEditor::timerCallback()
{
#ifdef CUSTOM_CMD
	if (needUpdate)
	{
		log.clear();
		int i = circularStrCurrent;
		int hist_num = 1;
		do
		{
			log.append(history[i], strlen(history[i]));
			hist_num++;
			i = (i + 1) % HISTORY_COUNT;
		} while (i != circularStrCurrent);
		// Write text to Juce
		consoleOutput->setText(log);
		// Scroll to the end
		consoleOutput->moveCaretToEnd(false);
		needUpdate = 0;
	}
#endif
}
void LiveProgEditor::paint(Graphics& g)
{
	g.fillAll(getLookAndFeel().findColour(ResizableWindow::backgroundColourId));
}
// ========== Define resized Method ==========
void LiveProgEditor::resized()
{
	const int width = getWidth();
	const int height = getHeight();
	fileComp->setBounds(5, 5, width - 10, 20);
	textContent->setBounds(5, 30, width - 185, height / 2 + 10);
	textContent->setColour(TextEditor::backgroundColourId, Colour(255, 255, 240));
	textContent->setColour(TextEditor::textColourId, Colour(0, 0, 0));
	table.setBounds(width - 185, 22, 187, height / 2 + 26);
	playButton.setBounds(5, height / 2 + (20 << 1) + 4, width / 8, 20);
//	errDisp->setBounds(width / 8 + 10, height / 2 + (20 << 1) + 4, width - width / 8 + 10, 20);

	juce::Rectangle<int> beforejSocial = getLocalBounds().reduced(1, 1).removeFromBottom(35);
	int asf = beforejSocial.getY() - height / 2 + (20 << 1) + 27 - 132;
	consoleOutput->setBounds(5, height / 2 + (20 << 1) + 27, width - 10, asf);
	consoleOutput->setColour(TextEditor::backgroundColourId, Colour(255, 255, 240));
	consoleOutput->setColour(TextEditor::textColourId, Colour(0, 0, 0));

	//jsocialButtons.setBounds(beforejSocial);
}
void LiveProgEditor::filenameComponentChanged(FilenameComponent* fileComponentThatHasChanged)
{
	if (fileComponentThatHasChanged == fileComp)
		readFile(fileComp->getCurrentFile());
}
void LiveProgEditor::readFile(const File& fileToRead)
{
	if (!fileToRead.existsAsFile())
		return;  // file doesn't exist
	FileInputStream inputStream(fileToRead);
	if (!inputStream.openedOk())
		return;  // failed to open
	textContent->clear();
	textContent->setColour(TextEditor::backgroundColourId, Colour(255, 255, 240));
	textContent->setText(inputStream.readString());
	textContent->pathOfCurrentFile = fileToRead.getFileName();
	processor.restoredText = textContent->getText();
}
