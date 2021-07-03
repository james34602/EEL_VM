#pragma once
#include <JuceHeader.h>
#include "PluginProcessor.h"
class Pop_Texteditor : public TextEditor
{
public:
	void addPopupMenuItems(PopupMenu& m, const MouseEvent*e)
	{
		TextEditor::addPopupMenuItems(m, e);
		m.addSeparator();
		m.addItem(57774, TRANS("Save file"));
	}
	bool isValidName(const String& filename)
	{
		if (filename.isEmpty()) return false;
		if (filename.trim() != filename) return false;
		return filename.removeCharacters("\\/:*?|<>\"") == filename;
	}
	void performPopupMenuAction(const int menuItemID)
	{
		TextEditor::performPopupMenuAction(menuItemID);
		if (menuItemID == 57774)
		{
			File resourceFile;
			if (pathOfCurrentFile.isNotEmpty() && !isValidName(pathOfCurrentFile))
				resourceFile = File(pathOfCurrentFile);
			FileChooser chooser("Saving liveprog program...", resourceFile, "*.eel;*.lua");
#ifndef __ANDROID__ 
			if (chooser.browseForFileToSave(true))
				resourceFile = chooser.getResult();
#endif
			pathOfCurrentFile = resourceFile.getFullPathName();
			FileOutputStream output(resourceFile);
			if (!output.openedOk())
				DBG("FileOutputStream didn't open correctly ...");
			output.setPosition(0);
			output.truncate();
			output.writeText(this->getText(), false, false, "\r\n");
			output.flush();
			if (output.getStatus().failed())
				DBG("An error occurred in the FileOutputStream");
			else
				fileCompPtr->setCurrentFile(resourceFile, true, NotificationType::dontSendNotification);
		}
	}
	String pathOfCurrentFile;
	FilenameComponent *fileCompPtr;
};
#define MAXITEMS 4096
#define UPDATEINTERVAL 1000
class TableComponent : public Component, public TableListBoxModel, private Timer
{
public:
	String displayStringAry[2][MAXITEMS + 1];
	char   ***ptrVar_name;
	float **ptrVar_value;
	int *ptrNB;
	int items;
	TableComponent(char ***nm, float **va, int *pnb)
	{
		items = 0;
		ptrVar_name = nm;
		ptrVar_value = va;
		ptrNB = pnb;
		addAndMakeVisible(table);
		table.setColour(ListBox::outlineColourId, Colours::grey);
		table.setOutlineThickness(1);
		// Add data
		table.getHeader().addColumn("Name", 1, 60, 40, -1, 31, TableHeaderComponent::defaultFlags);
		table.getHeader().addColumn("Value", 2, 50, 35, -1, 31, TableHeaderComponent::defaultFlags);
		table.setMultipleSelectionEnabled(true);
		startTimer(UPDATEINTERVAL);
		allowUpdate = 1;
	}
	~TableComponent()
	{
		stopTimer();
	}
	int getNumRows() override
	{
		return items;
	}
	int allowUpdate;
	void timerCallback() override
	{
		if (allowUpdate)
		{
			int counter = 0;
			for (int i = 0; i < *ptrNB; i++)
			{
				for (int j = 0; j < NSEEL_VARS_PER_BLOCK; j++)
				{
					displayStringAry[0][counter] = ptrVar_name[i][j];
					displayStringAry[1][counter] = String(ptrVar_value[i][j]);
					counter++;
				}
			}
			items = counter >= MAXITEMS ? MAXITEMS : counter;
			table.updateContent();
		}
		else
		{
			if (getTimerInterval() == 5000)
			{
				allowUpdate = 1;
				startTimer(UPDATEINTERVAL);
			}
			else
				startTimer(5000);
		}
	}
	void paintRowBackground(Graphics& g, int rowNumber, int /*width*/, int /*height*/, bool rowIsSelected) override
	{
		auto alternateColour = getLookAndFeel().findColour(ListBox::backgroundColourId)
			.interpolatedWith(getLookAndFeel().findColour(ListBox::textColourId), 0.03f);
		if (rowIsSelected)
			g.fillAll(Colours::lightblue);
		else if (rowNumber % 2)
			g.fillAll(alternateColour);
	}
	void paintCell(Graphics& g, int rowNumber, int, int width, int height, bool rowIsSelected) override
	{
		g.setColour(rowIsSelected ? Colours::darkblue : getLookAndFeel().findColour(0x100f006));
		g.setFont(font);
		String text = displayStringAry[0][rowNumber];
		g.drawText(text, 2, 0, width - 4, height, Justification::centredLeft, true);
		g.setColour(getLookAndFeel().findColour(ListBox::backgroundColourId));
		g.fillRect(width - 1, 0, 1, height);
	}
	Component* refreshComponentForCell(int rowNumber, int columnId, bool, Component* existingComponentToUpdate) override
	{
		if (columnId == 2)
		{
			EditableTextCustomComponent *textLabel = static_cast<EditableTextCustomComponent*> (existingComponentToUpdate);
			if (textLabel == nullptr)
				textLabel = new EditableTextCustomComponent(*this);
			textLabel->setRowAndColumn(rowNumber, columnId);
			return textLabel;
		}
		jassert(existingComponentToUpdate == nullptr);
		return nullptr;
	}
	String getText(const int columnNumber, const int rowNumber) const
	{
		return displayStringAry[columnNumber - 1][rowNumber];
	}
	void setText(const int columnNumber, const int rowNumber, const String& newText)
	{
		displayStringAry[columnNumber - 1][rowNumber] = newText;
		int counter = 0;
		for (int i = 0; i < *ptrNB; i++)
		{
			for (int j = 0; j < NSEEL_VARS_PER_BLOCK; j++)
			{
				String strInLoop = displayStringAry[0][counter];
				String str2Chk = displayStringAry[columnNumber - 2][rowNumber];
				if (strInLoop == str2Chk)
					if (newText.containsOnly("0123456789.-"))
						ptrVar_value[i][j] = newText.getFloatValue();
				counter++;
			}
		}
	}
	//==============================================================================
	void resized() override
	{
		table.setBoundsInset(BorderSize<int>(8));
	}
private:
	TableListBox table{ {}, this };
	Font font{ 15.0f };
	//==============================================================================
	class EditableTextCustomComponent : public Label
	{
	public:
		EditableTextCustomComponent(TableComponent& td) : owner(td)
		{
			setEditable(false, true, false);
		}
		void mouseDown(const MouseEvent& event) override
		{
			owner.table.selectRowsBasedOnModifierKeys(row, event.mods, false);
			Label::mouseDown(event);
			if (event.getNumberOfClicks() == 2)
			{
				owner.allowUpdate = 0;
			}
		}
		void textWasEdited() override
		{
			owner.setText(columnId, row, getText());
		}
		void setRowAndColumn(const int newRow, const int newColumn)
		{
			row = newRow;
			columnId = newColumn;
			setText(owner.getText(columnId, row), dontSendNotification);
		}
	private:
		TableComponent& owner;
		int row, columnId;
	};
	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TableComponent)
};
class LiveProgEditor : public AudioProcessorEditor, public FilenameComponentListener, private Timer
{
public:
	// declare constractor and deconstractor
	LiveProgEditor(LiveProgrammableDSP &p);
	~LiveProgEditor();

	// declare default methods
	void paint(Graphics &) override;
	void resized() override;
	std::unique_ptr<Pop_Texteditor> textContent;
	std::unique_ptr<Pop_Texteditor> consoleOutput;
//	Label *errDisp;

private:
	// audio processor's pointer
	LiveProgrammableDSP &processor;

	void filenameComponentChanged(FilenameComponent* fileComponentThatHasChanged) override;
	void readFile(const File& fileToRead);
	FilenameComponent *fileComp;
	TextButton playButton;
	// sliders
	TableComponent table;
	//jamesSocialButtons jsocialButtons;
	SharedResourcePointer<TooltipWindow> jtooltipWindow;
	String log;
	void timerCallback() override;
	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(LiveProgEditor)
};