#!/usr/bin/osascript

on run argv

    set argList to {}
    repeat with arg in argv
        set end of argList to quoted form of arg & space
    end repeat

    tell application "Terminal"
        #activatei
        set tid to do script "/usr/bin/ssh " & argList as string
        set current settings of tid to settings set "Ocean"
        set bounds of window 1 to {0, 0, 600, 400}
        #repeat until application Terminal is not running
        #    delay 1
        #end repeat
        #close window 1
    end tell

end run
