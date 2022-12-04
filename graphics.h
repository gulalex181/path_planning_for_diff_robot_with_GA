#include <vector>
#include <map>

#include <matplot/matplot.h>

template<typename T>
void plot_2d(
    const std::vector<T>& x,
    const std::vector<T>& y,
    const std::map<std::string, std::string>& keywords
) {

    std::map<std::string, std::string> local_keywords = keywords;

    // for (const auto& [key, value] : local_keywords) {
    //     std::cout << '[' << key << "] = " << value << "; ";
    // }

    // std::cout << std::endl;

    /* Default values */

    std::string title = "";
    static std::vector<std::string> legend = {};
    std::string color = "";
    std::string linestyle = "";
    std::string markerstyle = "";
    std::string hold = "";

    /* Get user values */

    auto search = local_keywords.find("title");
    
    if (search != local_keywords.end()) {
        title = search->second;
        local_keywords.erase(search);
    }

    search = local_keywords.find("legend");
    
    if (search != local_keywords.end()) {
        legend.push_back(search->second);
        local_keywords.erase(search);
    }

    search = local_keywords.find("color");
    
    if (search != local_keywords.end()) {
        color = search->second;
        local_keywords.erase(search);
    }

    search = local_keywords.find("linestyle");
    
    if (search != local_keywords.end()) {
        
        // none,          // no line
        // solid_line,    // "-"
        // dashed_line,   // "--"
        // dotted_line,   // ":"
        // dash_dot_line, // "-."
        
        linestyle = search->second;
        local_keywords.erase(search);

    }

    search = local_keywords.find("markerstyle");
    
    if (search != local_keywords.end()) {

        // none,                       // "" -> gnuplot linetype -1
        // plus_sign,                  // "+" -> gnuplot linetype 1
        // circle,                     // "o" -> gnuplot linetype 6
        // asterisk,                   // "*" -> gnuplot linetype 3
        // point,                      // "." -> gnuplot linetype 7
        // cross,                      // "x" -> gnuplot linetype 2
        // square,                     // "s" / "square" -> gnuplot linetype 4 / 5
        // diamond,                    // "d" / "diamond" -> gnuplot linetype 12 / 13
        // upward_pointing_triangle,   // "^" -> gnuplot linetype 8 / 9
        // downward_pointing_triangle, // "v" -> gnuplot linetype 10 / 11
        // right_pointing_triangle,    // ">" -> gnuplot linetype (doest not exist)
        // left_pointing_triangle,     // "<" -> gnuplot linetype (does not exist)
        // pentagram,                  // "p" / "pentagram" -> gnuplot linetype 14 / 15
        // hexagram,                   // "h" / "hexagram" -> gnuplot linetype (does not exist)

        markerstyle = search->second;
        local_keywords.erase(search);

    }

    search = local_keywords.find("hold");
    
    if (search != local_keywords.end()) {
        hold = search->second;
        local_keywords.erase(search);
    }

    /* Plot */

    // Add plot title
    matplot::title(title);

    // Show legend
    matplot::legend(legend);
    
    // Show grid
    matplot::grid(matplot::on);

    // Plot graphic
    matplot::plot(x, y)->line_width(2)
        .color(color)
        .line_style(linestyle)
        .marker_style(markerstyle);

    if (hold == "on") {

        matplot::hold(matplot::on);

    } else if (hold == "off") {

        matplot::hold(matplot::off);

        // Show plot
        matplot::show();

        // Clear legend vector
        legend.clear();

    } else {

        // Clear legend vector
        legend.clear();

    }

}

template<typename T>
void plot_3d(
    const std::vector<T>& x,
    const std::vector<T>& y,
    const std::vector<T>& z,
    const std::map<std::string, std::string>& keywords
) {

    std::map<std::string, std::string> local_keywords = keywords;

    // for (const auto& [key, value] : local_keywords) {
    //     std::cout << '[' << key << "] = " << value << "; ";
    // }

    // std::cout << std::endl;

    /* Default values */

    std::string title = "";
    static std::vector<std::string> legend = {};
    std::string color = "";
    std::string linestyle = "";
    std::string markerstyle = "";
    std::string hold = "";

    /* Get user values */

    auto search = local_keywords.find("title");
    
    if (search != local_keywords.end()) {
        title = search->second;
        local_keywords.erase(search);
    }

    search = local_keywords.find("legend");
    
    if (search != local_keywords.end()) {
        legend.push_back(search->second);
        local_keywords.erase(search);
    }

    search = local_keywords.find("color");
    
    if (search != local_keywords.end()) {
        color = search->second;
        local_keywords.erase(search);
    }

    search = local_keywords.find("linestyle");
    
    if (search != local_keywords.end()) {
        
        // none,          // no line
        // solid_line,    // "-"
        // dashed_line,   // "--"
        // dotted_line,   // ":"
        // dash_dot_line, // "-."
        
        linestyle = search->second;
        local_keywords.erase(search);

    }

    search = local_keywords.find("markerstyle");
    
    if (search != local_keywords.end()) {

        // none,                       // "" -> gnuplot linetype -1
        // plus_sign,                  // "+" -> gnuplot linetype 1
        // circle,                     // "o" -> gnuplot linetype 6
        // asterisk,                   // "*" -> gnuplot linetype 3
        // point,                      // "." -> gnuplot linetype 7
        // cross,                      // "x" -> gnuplot linetype 2
        // square,                     // "s" / "square" -> gnuplot linetype 4 / 5
        // diamond,                    // "d" / "diamond" -> gnuplot linetype 12 / 13
        // upward_pointing_triangle,   // "^" -> gnuplot linetype 8 / 9
        // downward_pointing_triangle, // "v" -> gnuplot linetype 10 / 11
        // right_pointing_triangle,    // ">" -> gnuplot linetype (doest not exist)
        // left_pointing_triangle,     // "<" -> gnuplot linetype (does not exist)
        // pentagram,                  // "p" / "pentagram" -> gnuplot linetype 14 / 15
        // hexagram,                   // "h" / "hexagram" -> gnuplot linetype (does not exist)

        markerstyle = search->second;
        local_keywords.erase(search);

    }

    search = local_keywords.find("hold");
    
    if (search != local_keywords.end()) {
        hold = search->second;
        local_keywords.erase(search);
    }

    /* Plot */

    // Add plot title
    matplot::title(title);

    // Show legend
    matplot::legend(legend);
    
    // Show grid
    matplot::grid(matplot::on);

    // Plot graphic
    matplot::plot3(x, y, z)->line_width(2)
        .color(color)
        .line_style(linestyle)
        .marker_style(markerstyle);

    if (hold == "on") {

        matplot::hold(matplot::on);

    } else if (hold == "off") {

        matplot::hold(matplot::off);

        // Show plot
        matplot::show();

        // Clear legend vector
        legend.clear();

    }

}