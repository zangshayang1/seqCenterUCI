// JS functions

function adminLoginCheck() {
    var username = document.forms["adminLoginForm"]["adminUsername"].value;
    var password = document.forms["adminLoginForm"]["adminPassword"].value;
    if (username != "jennyWu" || password != "jennyWu") {
        alert("You are not authorized to access the database.");
        return false;
    }
    return true;
}
