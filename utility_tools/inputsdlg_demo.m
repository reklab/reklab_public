%INPUTSDLG DEMO (Enhanced input dialog box with multiple data types)

% Written by: Takeshi Ikuma
% Date: Nov. 16 2009

clear; close all;

dim = [15 4];

Title = 'INPUTSDLG Demo Dialog';

% static text
Prompt = cell(1,2);
Prompt{1,1} = ['This demo illustrates every type of control that can be placed by INPUTSDLG function ' ... 
   'and demonstrates how Formats input can be used to layout these controls.'];
Formats(1,1).type = 'text';
Formats(1,1).size = [-1 0];
for k = 2:dim(2) % span the static text across the entire dialog
   Formats(1,k).type = 'none';
   Formats(1,k).limits = [0 1]; % extend from left
end

Prompt(2,:) = {'Bidder''s Name', 'Name'};
Formats(2,1).type = 'edit';
Formats(2,1).format = 'text';
Formats(2,1).size = 200; % automatically assign the height

Prompt(3,:) = {'Bidder''s SSN (no space or hyphen)', 'SSN'};
Formats(2,2).type = 'edit';
Formats(2,2).format = 'integer';
Formats(2,2).limits = [0 999999999]; % 9-digits (positive #)
Formats(2,2).size = 60;

Prompt(4,:) = {'Bidding Price', 'Price'};
Formats(2,3).type = 'edit';
Formats(2,3).format = 'float';
Formats(2,3).limits = [0 inf]; % non-negative decimal number

Prompt(5,:) = {'Enable enhanced mode.' 'EnableEnhancedMode'};
Formats(2,4).type = 'check';

Prompt(6,:) = {'Bidder''s Bio File','BioFile'};
Formats(3,1).type = 'edit';
Formats(3,1).format = 'file';
Formats(3,1).items = {'*.bio','Biography File (*.bio)';'*.*','All Files'};
Formats(3,1).limits = [0 1]; % single file get
Formats(3,1).size = [-1 0];

for k = 2:dim(2)-1 % span the file edit
   Formats(3,k).type = 'none';
   Formats(3,k).limits = [0 1]; % extend from left
end

Prompt(7,:) = {'Action','Action'};
Formats(3,4).type = 'list';
Formats(3,4).style = 'togglebutton';
Formats(3,4).items = {'Bid';'Decline';'Pass'};
for k = 4:6 % span the file edit vertically
   Formats(k,4).type = 'none';
   Formats(k,4).limits = [1 0]; % extend from above
end

Prompt(8,:) = {'Bidder''s Data Folder','DataFolder'};
Formats(4,1).type = 'edit';
Formats(4,1).format = 'dir';
Formats(4,1).size = [-1 0];

for k = 2:dim(2)-1 % span the dir edit
   Formats(4,k).type = 'none';
   Formats(4,k).limits = [0 1]; % extend from left
end

Prompt(9,:) = {'Save Bidding History To','SaveFile'};
Formats(5,1).type = 'edit';
Formats(5,1).format = 'file';
Formats(5,1).limits = [1 0]; % use uiputfile
Formats(5,1).size = [-1 0];

for k = 2:dim(2)-1 % span the file edit
   Formats(5,k).type = 'none';
   Formats(5,k).limits = [0 1]; % extend from left
end

Prompt(10,:) = {'Select Item Files','ItemFiles'};
Formats(6,1).type = 'edit';
Formats(6,1).format = 'file';
Formats(6,1).limits = [0 5]; % multi-select files
Formats(6,1).size = [-1 -1];
Formats(6,1).items = {'*.itm','Auction Item File';'*.*','All Files'};

for k = 2:dim(2)-1 % span the file edit
   Formats(6,k).type = 'none';
   Formats(6,k).limits = [0 1]; % extend from left
end

Prompt(11,:) = {'Choose Currency','MoneyUnit'};
Formats(7,1).type = 'list';
Formats(7,1).style = 'radiobutton';
Formats(7,1).items = {'U.S. Dollar' 'Euro';'Japanese Yen' ''};

Prompt(12,:) = {'Memo','Memo'};
Formats(7,2).type = 'edit';
Formats(7,2).format = 'text';
Formats(7,2).limits = [0 20]; % default: show 20 lines
Formats(7,2).size = [-1 0];
for k = 3:dim(2) % span horizontally
   Formats(7,k).type = 'none';
   Formats(7,k).limits = [0 1]; % extend from left
end
for row = 8:10
   for k = 3:dim(2) % span Memo area vertically
      Formats(row,k).type = 'none';
      Formats(row,k).limits = [1 0]; % extend from above
   end
end

Prompt(13,:) = {'Item Color','Color'};
Formats(8,1).type = 'list';
Formats(8,1).style = 'popupmenu';
Formats(8,1).items = {'Black','White','Red','Blue','Green','Yellow','Orange'};

Prompt(14,:) = {'Bidding Rate (left-slow;right-fast)','BidRate'};
Formats(9,1).type = 'range';
Formats(9,1).limits = [0 1];

Prompt(15,:) = {'Choose Auction Sites:','Site'};
Formats(10,1).type = 'list';
Formats(10,1).style = 'listbox';
Formats(10,1).items = {'www.auction1.com','www.auction2.com','www.bidme.com','www.bestvalu.com'};
Formats(10,1).limits = [0 2]; % multi-select

%%%% SETTING DIALOG OPTIONS
Options.WindowStyle = 'modal';
Options.Resize = 'on';
Options.Interpreter = 'tex';
Options.ApplyButton = 'on';

%%%% SETTING DEFAULT STRUCT
DefAns.Name = 'John Smith';
DefAns.SSN = 123456789;
DefAns.Price = 99.99;
DefAns.EnableEnhancedMode = true;
d = dir;
files = strcat([pwd filesep],{d(~[d.isdir]).name});
DefAns.BioFile = files{1};
DefAns.Action = 2; % = 'Decline'
DefAns.DataFolder = pwd;
DefAns.SaveFile = files{2};
DefAns.ItemFiles = files(3:end);
DefAns.MoneyUnit = 3; % yen
DefAns.Memo = 'Lorem ipsum dolor sit amet, consectetur adipiscing elit. Quisque elementum, dui sed sagittis vulputate, nulla tellus dapibus velit, pretium molestie odio lorem eget nisl. Cras volutpat gravida neque, vitae dictum eros aliquet non. Nam in sem sapien, sed condimentum justo. Cum sociis natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus. Fusce ut dignissim elit. Aenean ac massa arcu. Pellentesque rhoncus fermentum tortor et auctor. Proin eu velit dolor, eu semper nunc. Praesent in lacinia orci. Etiam ac diam nibh. Ut aliquet sapien sed metus blandit vitae aliquam enim ullamcorper. In nec ligula id quam interdum porttitor. Sed tempus turpis ut est convallis et interdum lectus dapibus. In dictum pulvinar erat id ornare.';
DefAns.Color = 5;
DefAns.BidRate = 0.75;
DefAns.Site = [1 2 4];

[Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options)
