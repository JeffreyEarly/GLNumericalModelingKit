Steps to make a new application based on GLNumericalModelingKit
===============================================================

1. Create a new Xcode project (command-line tool) with a git repository 'on my Mac'. The new project should be created in folder that also contains a checked-out copy of GLNumericalModelingKit.

2. Create a new Xcode workspace, and save it in the same directory as the new Xcode project. This new workspace will be empty, so drag in the Xcode project file for the new application as well as the Xcode project file for GLNumericalModelingKit.

3. Click on GLNumericalModelingKit and change the Location from 'Relative to Group' to 'Relative to Workspace.' Do the same thing for the application project as well.

4. Build GLNumericalModelingKit. You must do this first before proceeding to the next step, otherwise Xcode will reference the incorrect path.

5. In the 'Build Phases' tab of the application project, click the plus (+) button under the 'Link Binary With Libraries' and add GLNumericalModelingKit.dylib (or .a, or.framework). The GLNumericalModelingKit.dylib library should now appear under the 'Frameworks' group of the application project. Highlight the library and change its location from 'Relative to Group' to 'Relative to Build Products'.

6. Add "$(BUILT_PRODUCTS_DIR)/usr/local/include" to the header search paths under 'Build Settings' of the application project.

7. The application should now build correctly.

8. Add the workspace file to the local repository, `git add workspacename.xcworkspace`, then commit the local changes, `git commit -m 'some message'`.

9. Upload those changes to github (or where ever), using the standard methods. `git remote add github https://github.com/username/Hello-World.git` followed by `git push github master`.
