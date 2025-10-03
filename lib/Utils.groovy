class Utils {
    static String sanitize(String name) {
        return name.replaceAll('[^A-Za-z0-9_-]', '_').toLowerCase()
    }
}